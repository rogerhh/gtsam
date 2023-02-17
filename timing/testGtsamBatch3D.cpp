#include <gtsam/base/timing.h>
#include <gtsam/slam/dataset.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/sam/BearingRangeFactor.h>
#include <gtsam/geometry/Pose2.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/ISAM2.h>
#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/base/timing.h>
#include <cassert>


#include <unistd.h>
#include <getopt.h>


using namespace std;
using namespace gtsam;
using namespace gtsam::symbol_shorthand;

typedef Pose3 Pose;
typedef Point3 Point;

typedef NoiseModelFactor1<Pose> NM1;
typedef NoiseModelFactor2<Pose,Pose> NM2;
typedef BearingRangeFactor<Pose,Point> BR;


double chi2_red(const gtsam::NonlinearFactorGraph& graph, const gtsam::Values& config) {
    // Compute degrees of freedom (observations - variables)
    // In ocaml, +1 was added to the observations to account for the prior, but
    // the factor graph already includes a factor for the prior/equality constraint.
    //  double dof = graph.size() - config.size();
    int graph_dim = 0;
    for(const std::shared_ptr<gtsam::NonlinearFactor>& nlf: graph) {
        graph_dim += nlf->dim();
    }
    double dof = graph_dim - config.dim(); // kaess: changed to dim
    return 2. * graph.error(config) / dof; // kaess: added factor 2, graph.error returns half of actual error
}

int main(int argc, char *argv[]) {
    string dataset_name;
    int K = 1;
    double epsilon = 0.01;
    double d_error = 0.001;
    int max_iter = 10;

    // Get experiment setup
    static struct option long_options[] = {
        {"dataset", required_argument, 0, 'f'},
        {"K", required_argument, 0, 'k'},
        {"epsilon", required_argument, 0, 'e'},
        {"max_iter", required_argument, 0, 'm'},
        {"d_error", required_argument, 0, 'd'},
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
            default:
                cerr << "Unrecognized option" << endl;
                exit(1);
        }
    }

    assert(K > 0);

    cout << "Batch optimization" << endl
         << "Loading " << dataset_name << endl
         << "Parameters: K = " << K << ", epsilon = " << epsilon 
         << ", max_optimization_iter = " << max_iter << ", opt_stop_cond = " << d_error << endl;

    string datasetFile = findExampleDataFile(dataset_name);
    std::pair<NonlinearFactorGraph::shared_ptr, Values::shared_ptr> data =
        readG2o(datasetFile, true);

    NonlinearFactorGraph measurements = *data.first;
    Values initial = *data.second;
    Values estimate;
    NonlinearFactorGraph nfg;
    vector<int> opt_times;
    
    cout << "measurements.size() = " << measurements.size() << endl;
    cout << "values.size() = " << initial.size() << " ";

    cout << "Playing forward time steps..." << endl;

    try {
        int K_count = 0;
        int nextMeasurement = 0;
        for(size_t step=1; nextMeasurement < measurements.size(); ++step) {

            // Collect measurements and new variables for the current step
            if(step == 1) {
                //      cout << "Initializing " << 0 << endl;
                estimate.insert(0, Pose());
                // Add prior
                nfg.addPrior(0, Pose(), noiseModel::Unit::Create(6));
            }
            while(nextMeasurement < measurements.size()) {
                NonlinearFactor::shared_ptr measurementf = measurements[nextMeasurement];

                if(BetweenFactor<Pose>::shared_ptr measurement =
                        std::dynamic_pointer_cast<BetweenFactor<Pose> >(measurementf))
                {

                    // Stop collecting measurements that are for future steps
                    if(measurement->key1() > step || measurement->key2() > step)
                        break;

                    // Require that one of the nodes is the current one
                    if(measurement->key1() != step && measurement->key2() != step) {
                        // cout << measurement->key1() << " " << measurement->key2() << endl;

                        throw runtime_error("Problem in data file, out-of-sequence measurements");
                    }

                    // Add a new factor
                    nfg.push_back(measurement);

                    // Initialize the new variable
                    if(measurement->key1() == step && measurement->key2() == step - 1) {
                        if(step == 1) {
                            estimate.insert(step, measurement->measured().inverse());
                        }
                        else {
                            // Pose prevPose = isam2.calculateEstimate<Pose>(step-1);
                            assert(estimate.size() == step);
                            Pose prevPose = estimate.at<Pose>(step - 1);
                            estimate.insert(step, prevPose * measurement->measured().inverse());
                        }
                    } else if(measurement->key2() == step && measurement->key1() == step - 1) {
                        if(step == 1) {
                            estimate.insert(step, measurement->measured());
                        }
                        else {
                            // Pose prevPose = isam2.calculateEstimate<Pose>(step-1);
                            assert(estimate.size() == step);
                            Pose prevPose = estimate.at<Pose>(step - 1);
                            estimate.insert(step, prevPose * measurement->measured());
                        }
                    }

                }
                else
                {
                    throw std::runtime_error("Unknown factor type read from data file");
                }
                ++ nextMeasurement;
            }

            // if loop closure, run optimizer
            int iter = 0;
            int d = 0;
            if(K_count >= K || nextMeasurement == measurements.size()) {

                // LevenbergMarquardtOptimizer optimizer(nfg, estimate);
                // auto ordering = Ordering::Create(Ordering::NATURAL, nfg);
                auto ordering = Ordering::Create(Ordering::COLAMD, nfg);
                LevenbergMarquardtOptimizer optimizer(nfg, estimate, ordering);
                double lastError = optimizer.error();

                do {
                    auto opt_start = chrono::high_resolution_clock::now();
                    optimizer.iterate();
                    auto opt_end = chrono::high_resolution_clock::now();
                    d += chrono::duration_cast<chrono::microseconds>(opt_end - opt_start).count();
                    // cout << "Error: " << optimizer.error() << ", lambda: " << optimizer.lambda() << endl;
                    // break;    // Roger: This is weird. Why is this here?
                    double chi2 = chi2_red(nfg, estimate);
                    estimate = optimizer.values();
                    double curError = optimizer.error();
                    cout << "step = " << step << ", Chi2 = " << chi2 << ", eps = " << epsilon 
                         << ", graph_error = " << curError << endl;
                    if(chi2 <= epsilon) {
                        break;
                    }
                    iter++;
                    if(abs(lastError - curError) <= d_error) {
                        break;
                    }
                    else if(iter >= max_iter) {
                        cout << "Nonlinear optimization exceed max iterations: " 
                             << iter << " >= " << max_iter << ", chi2 = " << chi2 << endl;
                        break;
                    }
                    lastError = curError;
                } while(true);
                /*} while(!checkConvergence(optimizer.params().relativeErrorTol,
                            optimizer.params().absoluteErrorTol, optimizer.params().errorTol,
                            lastError, optimizer.error(), optimizer.params().verbosity));*/

                K_count = 0;
            }
            opt_times.push_back(d);
            K_count++;

        }

        double chi2 = chi2_red(nfg, estimate);
        cout << "final_chi2 = " << chi2 
             << ", final_error = " << nfg.error(estimate) << endl;

        for(int i = 0; i < opt_times.size(); i++) {
            cout << "step = " << i << ", opt_time = " << opt_times[i] << " us" << endl;
        }

        estimate.print();

    } catch(std::exception& e) {
        cout << e.what() << endl;
        return 1;
    }

    tictoc_print2_();

    return 0;
}
