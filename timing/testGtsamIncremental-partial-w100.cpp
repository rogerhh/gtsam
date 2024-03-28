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
    int relinearize_skip = 1;
    int print_frequency = 100;
    double epsilon = 0.01;
    double d_error = 0.001;
    int max_iter = 10;
    int num_steps = 1000000;
    int vio_tree_depth = INT_MAX;
    int relin_tree_depth = INT_MAX;
    int max_relin_factors = INT_MAX;
    int vio_relin_keys = INT_MAX;
    int max_relin_keys = INT_MAX;
    double relin_thresh = 0.1;
    int num_vio_step = 20;

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
        {"vio_tree_depth", required_argument, 0, 51},
        {"relin_tree_depth", required_argument, 0, 52},
        {"max_relin_factors", required_argument, 0, 53},
        {"vio_relin_keys", required_argument, 0, 54},
        {"max_relin_keys", required_argument, 0, 55},
        {"relin_thresh", required_argument, 0, 56},
        {"num_vio_step", required_argument, 0, 57},
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
            case 51:
                vio_tree_depth = atoi(optarg);
                break;
            case 52:
                relin_tree_depth = atoi(optarg);
                break;
            case 53:
                max_relin_factors = atoi(optarg);
                break;
            case 54:
                vio_relin_keys = atoi(optarg);
                break;
            case 55:
                max_relin_keys = atoi(optarg);
                break;
            case 56:
                relin_thresh = atof(optarg);
                break;
            case 57:
                num_vio_step = atoi(optarg);
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
         << ", relin_thresh = " << relin_thresh 
         << ", relin_tree_depth = " << relin_tree_depth 
         << ", vio_relin_keys = " << vio_relin_keys 
         << ", max_relin_keys = " << max_relin_keys 
         << ", num_vio_step = " << num_vio_step 
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

    ISAM2Params isam2params;
    isam2params.relinearizeThreshold = relin_thresh;
    isam2params.relinearizeSkip = relinearize_skip;
    ISAM2 isam2(isam2params);

    vector<int> update_times, calc_times;

    size_t nextMeasurement = 0;
    int K_count = 0;
    int print_count = 0;
    double last_chi2 = 0.0;
    int num_relin_keys = vio_relin_keys;
    Pose prevPose;
    Values newVariables;
    NonlinearFactorGraph newFactors;
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
            else {
                throw std::runtime_error("Unknown factor type read from data file");
            }
            ++ nextMeasurement;
        }

        // Update iSAM2
        int d1 = 0, d2 = 0;
        if(K_count == K || nextMeasurement == measurements.size()) {
            if(step % print_frequency == 0) {
                cout << "step = " << step << endl;
            }

            ISAM2UpdateParams params;
            params.force_relinearize = true;
            K_count = 0;
            Values estimate;
            auto start = chrono::high_resolution_clock::now();
            isam2.update(newFactors, newVariables, params, vio_tree_depth, vio_tree_depth, num_relin_keys);
            auto update_end = chrono::high_resolution_clock::now();
            auto calc_end = chrono::high_resolution_clock::now();
            d1 += chrono::duration_cast<chrono::microseconds>(update_end - start).count();
            d2 += chrono::duration_cast<chrono::microseconds>(calc_end - update_end).count();

            if(step >= num_vio_step) {
              num_relin_keys = max_relin_keys;
              estimate = isam2.calculateEstimate();
              last_chi2 = chi2_red(isam2.getFactorsUnsafe(), estimate);
              cout << "step = " << step << ", Chi2 = " << last_chi2 
                << ", graph_error = " << isam2.getFactorsUnsafe().error(estimate) << endl;
              estimate.print();
            }


            if(step >= num_steps) {
              break;
            }


            newVariables.clear();
            newFactors = NonlinearFactorGraph();
        }
        K_count++;
        update_times.push_back(d1);
        calc_times.push_back(d2);

    }

    Values estimate(isam2.calculateEstimate());
    double chi2 = chi2_red(isam2.getFactorsUnsafe(), estimate);
    cout << "final_chi2 = " << chi2 
         << ", final_error = " << isam2.getFactorsUnsafe().error(estimate) << endl;

    for(int i = 0; i < update_times.size(); i++) {
        cout << "step = " << i << ", update_time = " << update_times[i] << " us"
             << ", calc_time = " << calc_times[i] << " us" << endl;
    }

    tictoc_print2_();

}
