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
    int relinearize_skip = 4;
    int print_frequency = 100;
    double epsilon = 0.01;
    double d_error = 0.001;
    int max_iter = 10;
    int num_steps = 1000000;
    double relin_thresh = 0.1;
    string num_threads_infile = "";
    int num_threads = 1;
    string dataset_outdir = "";
    bool print_dataset = false;
    bool print_pred = false;
    bool print_traj = false;
    bool print_values = false;
    NoiseFormat noiseFormat = gtsam::NoiseFormatAUTO;

    // Get experiment setup
    static struct option long_options[] = {
        {"dataset", required_argument, 0, 'f'},
        {"K", required_argument, 0, 'k'},
        {"epsilon", required_argument, 0, 'e'},
        {"max_iter", required_argument, 0, 'm'},
        {"d_error", required_argument, 0, 'd'},
        {"relin_thresh", required_argument, 0, 'r'},
        {"relinearize_skip", required_argument, 0, 's'},
        {"print_frequency", required_argument, 0, 'p'},
        {"num_steps", required_argument, 0, 't'},
        {"num_threads", required_argument, 0, 50},
        {"dataset_outdir", required_argument, 0, 51},
        {"print_values", no_argument, 0, 52},
        {"num_threads_infile", required_argument, 0, 53},
        {"print_dataset", no_argument, 0, 150},
        {"print_pred", no_argument, 0, 151},
        {"print_traj", no_argument, 0, 152},
        {"noise_format", required_argument, 0, 54},
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
            case 'r':
                relin_thresh = atof(optarg);
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
            case 50:
                num_threads = atoi(optarg);
                break;
            case 51:
                dataset_outdir = string(optarg);
                break;
            case 52:
                print_values = true;
                break;
            case 53:
                num_threads_infile = string(optarg);
                break;
            case 54: {
                string noiseFormatString = string(optarg);
                if(noiseFormatString == "g2o") {
                  noiseFormat = NoiseFormatG2O;
                }
                else if(noiseFormatString == "toro") {
                  noiseFormat = NoiseFormatTORO;
                }
                else if(noiseFormatString == "graph") {
                  noiseFormat = NoiseFormatGRAPH;
                }
                else if(noiseFormatString == "cov") {
                  noiseFormat = NoiseFormatCOV;
                }
                else if(noiseFormatString == "auto") {
                  noiseFormat = NoiseFormatAUTO;
                }
                else {
                  noiseFormat = NoiseFormatAUTO;
                }
                break;
            }
            case 150:
                print_dataset = true;
                break;
            case 151:
                print_pred = true;
                break;
            case 152:
                print_traj = true;
                break;
            default:
                cerr << "Unrecognized option" << endl;
                exit(1);
        }
    }

    bool use_num_threads_file = false;
    ifstream num_threads_fin;
    if(num_threads_infile != "") {
        use_num_threads_file = true;
        num_threads_fin.open(num_threads_infile);
        if(!num_threads_fin.is_open()) {
            cout << "Error opening file: " << num_threads_infile << endl;
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
        load2D(datasetFile, SharedNoiseModel(), 0, false, true, noiseFormat);

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

            if(use_num_threads_file) {
                num_threads_fin >> num_threads;
            }

            isam2.update_resource_aware(newFactors, newVariables, params, num_threads);
            auto update_end = chrono::high_resolution_clock::now();
            estimate = isam2.calculateEstimate();
            auto calc_end = chrono::high_resolution_clock::now();
            d1 += chrono::duration_cast<chrono::microseconds>(update_end - start).count();
            d2 += chrono::duration_cast<chrono::microseconds>(calc_end - update_end).count();

            if(step % print_frequency == 0) {
                // estimate = isam2.calculateEstimate();
                // cout << "Theta = " << endl;
                // estimate.print();
            }
            if(step >= num_steps) {
                break;
            }

            last_chi2 = chi2_red(isam2.getFactorsUnsafe(), estimate);
            cout << "chi2 = " << last_chi2 << endl;

            // last_chi2 = chi2_red(isam2.getFactorsUnsafe(), estimate);
            // print_count++;
            // if(print_frequency != 0 && print_count % print_frequency == 0) {
            //     cout << "step = " << step << ", Chi2 = " << last_chi2 
            //          << ", graph_error = " << isam2.getFactorsUnsafe().error(estimate) << endl;
            // }

            // if(K > 1) {
            //     NonlinearFactorGraph dummy_nfg;
            //     Values dummy_vals;
            //     int iter = 0;

            //     while(1) {
            //         if(last_chi2 <= epsilon) {
            //             break;
            //         }
            //         auto start = chrono::high_resolution_clock::now();
            //         isam2.update(dummy_nfg, dummy_vals);
            //         auto update_end = chrono::high_resolution_clock::now();
            //         estimate = isam2.calculateEstimate();
            //         auto calc_end = chrono::high_resolution_clock::now();
            //         d1 += chrono::duration_cast<chrono::microseconds>
            //                     (update_end - start).count();
            //         d2 += chrono::duration_cast<chrono::microseconds>
            //                     (calc_end - update_end).count();

            //         double chi2 = chi2_red(isam2.getFactorsUnsafe(), estimate);

            //         if(print_frequency != 0 && print_count % print_frequency == 0) {
            //             cout << "step = " << step << ", Chi2 = " << last_chi2 
            //                 << ", graph_error = " << isam2.getFactorsUnsafe().error(estimate) << endl;
            //         }

            //         if(abs(last_chi2 - chi2) < d_error) {
            //             break;
            //         }

            //         last_chi2 = chi2;
            //         iter++;
            //         if(iter >= max_iter) {
            //             cout << "Nonlinear optimization exceed max iterations: " 
            //                  << iter << " >= " << max_iter << ", chi2 = " << chi2 << endl;
            //             break;
            //         }
            //     }
            // }
            newVariables.clear();
            newFactors = NonlinearFactorGraph();

            if(step % print_frequency == 0) {
                if(dataset_outdir != "") {
                    if(print_dataset) {
                        string outfile = dataset_outdir + "/step-" + to_string(step) + ".out";
                        ofstream fout(outfile);

                        if(!fout.is_open()) {
                            cerr << "Cannot open file: " << outfile << endl;
                            exit(1);
                        }

                        isam2.extractFullTree(fout, print_values);
                        if(print_values) {
                            isam2.extractDelta(fout);
                        }
                    }

                    if(print_pred) {
                        string pred_outfile = dataset_outdir + "/step-" + to_string(step) + "-pred_cycles.out";
                        ofstream pred_fout(pred_outfile);

                        if(!pred_fout.is_open()) {
                            cerr << "Cannot open file: " << pred_outfile << endl;
                            exit(1);
                        }

                        int num_threads = 1;
                        isam2.extractPredictedCycles(pred_fout, num_threads);
                    }

                    if(print_traj) {
                        string traj_outfile = dataset_outdir + "/step-" + to_string(step) + "_traj.txt";

                        ofstream traj_fout(traj_outfile);

                        if(!traj_fout.is_open()) {
                            cerr << "Cannot open file: " << traj_outfile << endl;
                            exit(1);
                        }

                        estimate.print_kitti_pose2(traj_fout);
                    }
                }
            }
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

    estimate.print();

    tictoc_print2_();

}
