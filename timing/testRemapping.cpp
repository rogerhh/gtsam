#include <gtsam/slam/dataset.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/sam/BearingRangeFactor.h>
#include <gtsam/geometry/Pose2.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/ISAM2.h>
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
#include <cstdlib>

#include <unordered_set>

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

void recursiveMarkMarginalizeAffectedKeys(const Key& key, std::set<Key>& additionKeys) {

}

void recursiveMarkAffectedKeys(const Key& key,
    const ISAM2Clique::shared_ptr& clique, std::set<Key>& additionalKeys) {

  // Check if the separator keys of the current clique contain the specified key
  if (std::find(clique->conditional()->beginParents(),
      clique->conditional()->endParents(), key)
      != clique->conditional()->endParents()) {

    // Mark the frontal keys of the current clique
    for(Key i: clique->conditional()->frontals()) {
      additionalKeys.insert(i);
    }

    // Recursively mark all of the children
    for(const ISAM2Clique::shared_ptr& child: clique->children) {
      recursiveMarkAffectedKeys(key, child, additionalKeys);
    }
  }
  // If the key was not found in the separator/parents, then none of its children can have it either
}

void createOrderingConstraints(
    Key max_key,
    const KeyVector& marginalizableKeys,
    const set<Key>& activeKeys,
    FastMap<Key, int>& constrainedKeys) {
  if (marginalizableKeys.size() > 0) {
    constrainedKeys = FastMap<Key, int>();
    // Generate ordering constraints so that the marginalizable variables will be eliminated first
    // Set all variables to Group1
    for(Key k : activeKeys) {
      constrainedKeys[k] = 1;
    }
    // Set marginalizable variables to Group0
    for(Key key: marginalizableKeys) {
      constrainedKeys[key] = 0;
    }
  }
}


void isam2_marginalization(
    Key max_key, 
    const KeyVector& marginalizableKeys,
    const set<Key>& activeKeys,
    ISAM2& isam2,
    FastMap<Key, int>* constrainedKeys,
    KeyList* additionalMarkedKeys) {

  createOrderingConstraints(max_key, marginalizableKeys, activeKeys, *constrainedKeys);

  // Mark additional keys between the marginalized keys and the leaves
  std::set<Key> additionalKeys;
  for(Key key: marginalizableKeys) {
    ISAM2Clique::shared_ptr clique = isam2[key];
    for(const ISAM2Clique::shared_ptr& child: clique->children) {
      recursiveMarkAffectedKeys(key, child, additionalKeys);
    }
  }
  *additionalMarkedKeys = KeyList(activeKeys.begin(), activeKeys.end());
}

int main(int argc, char *argv[]) {

    string dataset_name;
    string marginalize_test_inputs;
    int K = 1;
    int relinearize_skip = 4;
    int print_frequency = 100;
    double epsilon = 0.01;
    double d_error = 0.001;
    int max_iter = 10;
    int num_steps = 1000000;
    bool is_new = true;
    int seed = 0;

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
        {"marginalize_test_inputs", required_argument, 0, 1},
        {"new", no_argument, 0, 2},
        {"old", no_argument, 0, 3},
        {"seed", required_argument, 0, 4},
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
            case 1:
                marginalize_test_inputs = string(optarg);
                break;
            case 2:
                is_new = true;
                break;
            case 3:
                is_new = false;
                break;
            case 4:
                seed = atoi(optarg);
                srand(seed);
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

    cout << "Marginalize test file: " << marginalize_test_inputs << endl;

    if(is_new) {
        cout << "Testing NEW code" << endl;
    }
    else {
        cout << "Testing OLD code" << endl;
    }

    ifstream margin_fin(marginalize_test_inputs);
    if(!margin_fin.is_open()) {
      cerr << "Error opening file: " << marginalize_test_inputs;
      exit(1);
    }

    // Read in marginalization test
    string op;
    int op_count = 0;
    int remove_or_marginalize_index = 0;
    vector<int> remove_marginalize_steps;
    vector<vector<int>> remove_factor_indices;
    vector<vector<int>> marginalize_keys;
    unordered_set<int> marginalized_keys; // keys already marginalized
    // Because of marginalization, some factors may not be added. 
    // We need a map to keep track of factors and actual factor indices
    unordered_map<int, int> real_factor_indices; 
    do {
      margin_fin >> op;
    } while(op != "START");
    while(margin_fin >> op >> op_count) {
      if(op != "STEP") {
        throw runtime_error("op is not STEP: " + op);
      } 
      remove_marginalize_steps.push_back(op_count);
      remove_factor_indices.push_back(vector<int>());
      marginalize_keys.push_back(vector<int>());
      margin_fin >> op >> op_count;
      if(op != "REMOVE") {
        throw runtime_error("op is not REMOVE: " + op);
      } 
      int k;
      for(int i = 0; i < op_count; i++) {
        margin_fin >> k;
        remove_factor_indices.back().push_back(k);
      }
      margin_fin >> op >> op_count;
      if(op != "MARGINALIZE") {
        throw runtime_error("op is not REMOVE: " + op);
      } 
      for(int i = 0; i < op_count; i++) {
        margin_fin >> k;
        marginalize_keys.back().push_back(k);
      }
    }

    string datasetFile = findExampleDataFile(dataset_name);
    // string datasetFile = findExampleDataFile("victoria_park");
    std::pair<NonlinearFactorGraph::shared_ptr, Values::shared_ptr> data =
        load2D(datasetFile);

    NonlinearFactorGraph measurements = *data.first;
    Values initial = *data.second;

    cout << "measurements.size() = " << measurements.size() << endl;
    cout << "values.size() = " << initial.size() << " ";

    cout << "Playing forward time steps..." << endl;

    ISAM2Params isam2params;
    isam2params.relinearizeSkip = relinearize_skip;
    ISAM2 isam2(isam2params);

    vector<int> update_times, calc_times;

    size_t nextMeasurement = 0;
    int K_count = 0;
    int print_count = 0;
    double last_chi2 = 0.0;
    Pose prevPose;
    Values newVariables;
    NonlinearFactorGraph newFactors;

    set<RemappedKey> activeKeys;
    activeKeys.insert(0);

    map<Key, RemappedKey> keyTransformMap;
    map<Key, RemappedKey> keyUntransformMap;

    for(size_t step=1; nextMeasurement < measurements.size(); ++step) {
      
      RemappedKey var;
      bool inserted;
      do {
        var = rand() % 100000;
        inserted = keyUntransformMap.insert({var, step}).second;
      } while(!inserted);

      keyTransformMap.insert({step, var});

        activeKeys.insert(var);

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
                auto key1 = measurement->key1();
                auto key2 = measurement->key2();
                auto min_key = min(key1, key2);

                if(marginalized_keys.find(min_key) != marginalized_keys.end()) {
                   cout << "Skipped factor: " << nextMeasurement << endl;
                }
                else {
                  int real_factor_index;
                  if(is_new) {
                    real_factor_index = isam2.getCholeskyEliminationTree().numFactors() + newFactors.size();
                  }
                  else {
                    real_factor_index = isam2.getFactorsUnsafe().size() + newFactors.size();
                  }
                  real_factor_indices.insert({nextMeasurement, real_factor_index});

                  newFactors.push_back(measurement);

                  // Initialize the new variable
                  if(measurement->key1() == step && measurement->key2() == step-1) {
                    Pose newPose;
                    if(step == 1) {
                      newPose = measurement->measured().inverse();
                    }
                    else {
                      if(isam2.valueExists(step - 1)) {
                        prevPose = isam2.calculateEstimate<Pose>(step - 1);
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
                      if(isam2.valueExists(step - 1)) {
                        prevPose = isam2.calculateEstimate<Pose>(step - 1);
                      }
                      newPose = prevPose * measurement->measured();
                    }
                    newVariables.insert(step, newPose);
                    prevPose = newPose;
                  }
                }
            }
            else {
                throw std::runtime_error("Unknown factor type read from data file");
            }
            ++ nextMeasurement;
        }

        // Update iSAM2
        if(step % print_frequency == 0) {
          cout << "step = " << step << endl;
        }

        KeyList additionalMarkedKeys;
        FactorIndices factorsToRemove;
        boost::optional<FastMap<Key, int> > constrainedKeys = FastMap<Key, int>();

        FastList<Key> leafKeys;

        if(!remove_marginalize_steps.empty() &&
            step == remove_marginalize_steps[remove_or_marginalize_index]) {
          for(auto factor_index : remove_factor_indices[remove_or_marginalize_index]) {
            auto real_factor_index = real_factor_indices[factor_index];
            factorsToRemove.push_back(real_factor_index);
            if(is_new) {
              assert(isam2.getCholeskyEliminationTree().nonlinearFactorAt(real_factor_index) == measurements[factor_index]);
            }
            else {
              assert(isam2.getFactorsUnsafe()[real_factor_index] == measurements[factor_index]);
            }
          }

          auto& marginalizableKeys = marginalize_keys[remove_or_marginalize_index];

          if(is_new) {
            std::set<Key> additionalKeys;

            if(!marginalizableKeys.empty()) {

              for(Key key : activeKeys) {
                constrainedKeys->operator[](key) = 1;
              }

              for(Key key : marginalizableKeys) {
                // Recursively mark children of marginalized keys
                isam2.getCholeskyEliminationTree().getAffectedDescendantKeys(key, additionalKeys);
                // Recursively mark ancestors of marginalized keys
                isam2.getCholeskyEliminationTree().getAffectedKeys(key, additionalKeys);
                constrainedKeys->operator[](key) = 0;
                leafKeys.push_back(key);
              }

            }

            additionalMarkedKeys = KeyList(additionalKeys.begin(), additionalKeys.end());
          }
          else {
            KeyVector marginalizableKeyVector(marginalizableKeys.begin(), marginalizableKeys.end());
            FastMap<Key, int> constrainedKeysMap;
            isam2_marginalization(step, marginalizableKeyVector, activeKeys, isam2, &constrainedKeysMap, &additionalMarkedKeys);
            constrainedKeys = constrainedKeysMap;
            assert(constrainedKeys.is_initialized());
            leafKeys = FastList<Key>(marginalizableKeys.begin(), marginalizableKeys.end());
          }

          remove_or_marginalize_index++;
        }

        if(constrainedKeys->empty()) {
          constrainedKeys = boost::none;
        }

        // Update iSAM2
        auto isamResult_ = isam2.update(newFactors, newVariables,
                        factorsToRemove, constrainedKeys, boost::none, additionalMarkedKeys, true);

        if(!leafKeys.empty()) {
          for(Key key : leafKeys) {
            activeKeys.erase(key);
            marginalized_keys.insert(key);
          }
          isam2.marginalizeLeaves(leafKeys, boost::none, boost::none);
        }


        if(step % print_frequency == 0) {
          auto estimate = isam2.calculateEstimate();
          cout << "Theta = " << endl;
          estimate.print();
        }
        if(step >= num_steps) {
          break;
        }

        newVariables.clear();
        newFactors = NonlinearFactorGraph();

    }

    Values estimate(isam2.calculateEstimate());
    estimate.print();

    isam2.getCholeskyEliminationTree().printOrderingUnmapped(cout);
    isam2.getCholeskyEliminationTree().printOrderingRemapped(cout);

}
