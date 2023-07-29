#include <gtsam/linear/CholeskyEliminationTreeFactorWrapper.h>
#include <gtsam/linear/CholeskyEliminationTreeNode.h>

#include <gtsam/slam/dataset.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/sam/BearingRangeFactor.h>
#include <gtsam/geometry/Pose2.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/base/timing.h>

#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <chrono>

#include <unistd.h>
#include <getopt.h>

using namespace std;
using namespace gtsam;
using namespace gtsam::symbol_shorthand;

typedef Pose2 Pose;

typedef NoiseModelFactor1<Pose> NM1;
typedef NoiseModelFactor2<Pose,Pose> NM2;
typedef BearingRangeFactor<Pose,Point2> BR;

double chi2_red(const gtsam::NonlinearFactorGraph& graph, const gtsam::Values& config) {
    // Compute degrees of freedom (observations - variables)
    // In ocaml, +1 was added to the observations to account for the prior, but
    // the factor graph already includes a factor for the prior/equality constraint.
    //  double dof = graph.size() - config.size();
    int graph_dim = 0;
    for(const boost::shared_ptr<gtsam::NonlinearFactor>& nlf: graph) {
        graph_dim += nlf->dim();
    }
    double dof = graph_dim - config.dim(); // kaess: changed to dim
    return 2. * graph.error(config) / dof; // kaess: added factor 2, graph.error returns half of actual error
}

int main(int argc, char *argv[]) {

    string dataset_name;
    int K = 1;
    int relinearize_skip = 4;
    int print_frequency = 100;
    double epsilon = 0.01;
    double d_error = 0.001;
    int max_iter = 10;
    int num_steps = 1000000;

    // Get experiment setup
    static struct option long_options[] = {
        {"dataset", required_argument, 0, 'f'},
        {"K", required_argument, 0, 'k'},
        {"epsilon", required_argument, 0, 'e'},
        {"max_iter", required_argument, 0, 'm'},
        {"d_error", required_argument, 0, 'd'},
        {"relinearize_skip", required_argument, 0, 's'},
        {"print_frequency", required_argument, 0, 'p'},
        {"num_steps", required_argument, 0, 't'},
        {0, 0, 0, 0}
    };
    int opt, option_index;
    while((opt = getopt_long(argc, argv, "f:k:e:m:d:", long_options, &option_index)) != -1) {
        switch(opt) {
            case 'f':
                dataset_name = string(optarg);
                break;
            case 'k':
                K = atoi(optarg);
                break;
            case 'e':
                epsilon = atof(optarg);
                break;
            case 'm':
                max_iter = atoi(optarg);
                break;
            case 'd':
                d_error = atof(optarg);
                break;
            case 's':
                relinearize_skip = atoi(optarg);
                break;
            case 'p':
                print_frequency = atoi(optarg);
                break;
            case 't':
                num_steps = atoi(optarg);
                break;
            default:
                cerr << "Unrecognized option" << endl;
                exit(1);
        }
    }

    cout << "Incremental optimization" << endl
         << "Loading " << dataset_name << endl
         << "Parameters: K = " << K << ", epsilon = " << epsilon 
         << ", max_optimization_iter = " << max_iter 
         << ", opt_stop_cond = " << d_error 
         << ", relinearize_skip = " << relinearize_skip 
         << ", print_frequency = " << print_frequency 
         << endl;

    string datasetFile = findExampleDataFile(dataset_name);
    // string datasetFile = findExampleDataFile("victoria_park");
    std::pair<NonlinearFactorGraph::shared_ptr, Values::shared_ptr> data =
        load2D(datasetFile);

    NonlinearFactorGraph measurements = *data.first;
    Values initial = *data.second;

    cout << "measurements.size() = " << measurements.size() << endl;
    cout << "values.size() = " << initial.size() << " ";

    cout << "Playing forward time steps..." << endl;

    vector<int> update_times, calc_times;

    CholeskyEliminationTree etree;
    size_t nextMeasurement = 0;
    int K_count = 0;
    int print_count = 0;
    double last_chi2 = 0.0;
    Pose prevPose;
    Values newVariables;
    NonlinearFactorGraph newFactors;
    NonlinearFactorGraph allFactors;

    for(size_t step=1; nextMeasurement < measurements.size(); ++step) {

        // Collect measurements and new variables for the current step
        if(step == 1) {
            newVariables.insert(0, Pose());
            // Add prior
            newFactors.addPrior(0, Pose(), noiseModel::Unit::Create(3));
        }
        while(nextMeasurement < measurements.size()) {

            NonlinearFactor::shared_ptr measurementf = measurements[nextMeasurement];

            if(BetweenFactor<Pose>::shared_ptr measurement =
                    boost::dynamic_pointer_cast<BetweenFactor<Pose> >(measurementf))
            {

                // Stop collecting measurements that are for future steps
                if(measurement->key1() > step || measurement->key2() > step)
                    break;

                // Require that one of the nodes is the current one
                if(measurement->key1() != step && measurement->key2() != step) {
                    cout << measurement->key1() << " " << measurement->key2() << endl;

                    throw runtime_error("Problem in data file, out-of-sequence measurements");
                }

                // Add a new factor
                newFactors.push_back(measurement);

                // Initialize the new variable
                if(measurement->key1() == step && measurement->key2() == step-1) {
                    Pose newPose;
                    if(step == 1) {
                        newPose = measurement->measured().inverse();
                    }
                    else {
                        if(step > 0) {
                          prevPose = newVariables.at<Pose>(step - 1);
                        }
                        newPose = prevPose * measurement->measured().inverse();
                    }
                    newVariables.insert(step, newPose);
                    prevPose = newPose;
                } else if(measurement->key2() == step && measurement->key1() == step-1) {
                    Pose newPose;
                    if(step == 1) {
                        newPose = measurement->measured();
                    }
                    else {
                        if(step > 0) {
                          prevPose = newVariables.at<Pose>(step - 1);
                        }
                        newPose = prevPose * measurement->measured();
                    }
                    newVariables.insert(step, newPose);
                    prevPose = newPose;
                }

            }
            else {
                throw std::runtime_error("Unknown factor type read from data file");
            }
            ++ nextMeasurement;
        }
        
        RemappedKeySet affectedKeys, observedKeys;
        FactorIndices newFactorIndices = allFactors.add_factors(newFactors, true);

        etree.addVariables(newVariables);
        etree.markAffectedKeys(newFactors, newFactorIndices, KeySet(), ISAM2UpdateParams(), &affectedKeys, &observedKeys);

        etree.testFactorWrapper(initial);

        newVariables.clear();
        newFactors = NonlinearFactorGraph();

        if(step >= num_steps) {
          break;;
        }
    }
}

void CholeskyEliminationTree::testFactorWrapper(const Values& theta) {
  /*
  BlockIndexVector blockIndices;
  vector<tuple<Key, size_t>> keyWidthPairs({{2, 3}, {0, 3}, {4, 3}, {1, 3}});
  size_t curCol = 0;
  for(auto&[key, width] : keyWidthPairs) {
    blockIndices.push_back({addRemapKey(key), curCol, width});
    curCol += width;
  }
  blockIndices.push_back({0, curCol, 1});

  BlockIndexMap blockMap;
  for(auto& p : blockIndices) {
    blockMap.insert({std::get<BLOCK_INDEX_KEY>(p), 
                    {get<BLOCK_INDEX_ROW>(p), get<BLOCK_INDEX_HEIGHT>(p)}});
  }

  size_t totalHeight = get<BLOCK_INDEX_ROW>(blockIndices.back()) + 1;
  Matrix m(totalHeight, totalHeight);
  m.setZero();

  struct Pred {
    CholeskyEliminationTree* etree = nullptr;
    Pred(CholeskyEliminationTree* etree_in) : etree(etree_in) {}
    bool operator()(const RemappedKey& key) const { 
      return key != 0; 
    }
  };

  for(sharedFactorWrapper factorWrapper : factors_) {
    factorWrapper->linearizeIfNeeded(theta);
    factorWrapper->getLinearFactor()->print();
    factorWrapper->updateHessian(m, -1, blockMap, Pred(this));
  }

  cout << "m = \n" << m << endl;
  */
  assert(0);
}

