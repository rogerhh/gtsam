#include "gtsam/linear/NoiseModel.h"
#include <gtsam/slam/dataset.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/sam/BearingRangeFactor.h>
#include <gtsam/geometry/Pose3.h>
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

typedef Pose3 Pose;

typedef NoiseModelFactor1<Pose> NM1;
typedef NoiseModelFactor2<Pose,Pose> NM2;
typedef BearingRangeFactor<Pose,Point3> BR;

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
    int vio_lag = 5;
    double relin_thresh = 0.1;
    NoiseFormat noiseFormat = gtsam::NoiseFormatAUTO;
    string dataset_outdir = "";
    bool run_lc = false;
    int lc_period = 30;
    int loop_size = lc_period;

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
        {"vio_lag", required_argument, 0, 52},
        {"dataset_outdir", required_argument, 0, 54},
        {"run_lc", no_argument, 0, 55},
        {"lc_period", required_argument, 0, 56},
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
            case 52:
                vio_lag = atoi(optarg);
                break;
            case 53: {
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
            case 54:
                dataset_outdir = string(optarg);
                break;
            case 55:
                run_lc = true;
                break;
            case 56: {
                lc_period = atoi(optarg);
                loop_size = lc_period;
                break;
            }
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
         << ", print_frequency = " << print_frequency 
         << endl;

    string datasetFile = findExampleDataFile(dataset_name);
    // string datasetFile = findExampleDataFile("victoria_park");
    std::pair<NonlinearFactorGraph::shared_ptr, Values::shared_ptr> data =
        readG2o(datasetFile, true);

    NonlinearFactorGraph measurements = *data.first;
    Values initial = *data.second;

    cout << "measurements.size() = " << measurements.size() << endl;
    cout << "values.size() = " << initial.size() << " ";

    cout << "Playing forward time steps..." << endl;


    // Keep a list of all factors up to that point
    // In the vio case, just use LC poses 

    ISAM2Params isam2params;
    isam2params.relinearizeThreshold = relin_thresh;
    isam2params.relinearizeSkip = relinearize_skip;
    shared_ptr<ISAM2> vio_isam2 = std::make_shared<ISAM2>(isam2params);
    ISAM2 isam2(isam2params);

    size_t nextMeasurement = 0;
    int K_count = 0;
    int print_count = 0;
    double last_chi2 = 0.0;
    Pose prevPose;
    Values newVariables;
    NonlinearFactorGraph newFactors;

    Pose vioPrevPose;
    Values vioNewVariables;
    NonlinearFactorGraph vioNewFactors;

    int loop_start_step = 0;
    bool new_loop = false;

    for(size_t step=1; nextMeasurement < measurements.size(); ++step) {

        // Collect measurements and new variables for the current step
        if(step == 1) {
            newVariables.insert(0, Pose());
            // Add prior
            newFactors.addPrior(0, Pose(), noiseModel::Unit::Create(6));

            vioNewVariables.insert(0, Pose());
            vioNewFactors.addPrior(0, Pose(), noiseModel::Unit::Create(6));
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

                Key k1;
                Key k2;
                
                if(measurement->key1() > measurement->key2()) {
                  k1 = measurement->key1();
                  k2 = measurement->key2();
                }
                else {
                  k1 = measurement->key2();
                  k2 = measurement->key1();
                }

                if(k1 - k2 <= vio_lag) {
                  vioNewFactors.push_back(measurement);
                }
                else {
                    cout << "discarding factor " << k1 << " " << k2 << " " << vio_lag << endl;
                }

                if(k1 - k2 >= loop_size) {
                    cout << "detected loop" << endl;
                    if(loop_start_step == 0) {
                        loop_start_step = step;
                    }
                    else {
                        cout << "new loop" << endl;
                        new_loop = true;
                    }
                }

                // Add a new factor
                newFactors.push_back(measurement);

                // Initialize the new variable
                if(!newVariables.exists(step)) {
                  if(measurement->key1() == step && measurement->key2() == step-1) {
                    Pose newPose;
                    Pose vioNewPose;
                    if(step == 1) {
                      newPose = measurement->measured().inverse();
                      vioNewPose = newPose;
                    }
                    else {
                      if(isam2.valueExists(step - 1)) {
                        prevPose = isam2.calculateEstimate<Pose>(step - 1);
                        vioPrevPose = vio_isam2->calculateEstimate<Pose>(step - 1);
                      }
                      newPose = prevPose * measurement->measured().inverse();
                      vioNewPose = vioPrevPose * measurement->measured().inverse();
                    }
                    newVariables.insert(step, newPose);
                    prevPose = newPose;

                    vioNewVariables.insert(step, vioNewPose);
                    vioPrevPose = vioNewPose;
                  } else if(measurement->key2() == step && measurement->key1() == step-1) {
                    Pose newPose;
                    Pose vioNewPose;
                    if(step == 1) {
                      newPose = measurement->measured();
                      vioNewPose = newPose;
                    }
                    else {
                      if(isam2.valueExists(step - 1)) {
                        prevPose = isam2.calculateEstimate<Pose>(step - 1);
                        vioPrevPose = vio_isam2->calculateEstimate<Pose>(step - 1);
                      }
                      newPose = prevPose * measurement->measured();
                      vioNewPose = vioPrevPose * measurement->measured();
                    }
                    newVariables.insert(step, newPose);
                    prevPose = newPose;

                    vioNewVariables.insert(step, vioNewPose);
                    vioPrevPose = vioNewPose;
                  }
                }

            }
            else {
                throw std::runtime_error("Unknown factor type read from data file");
            }
            ++ nextMeasurement;
        }

        bool lc_running = false;
        cout << "loop_start_step = " << loop_start_step << " " << step << endl;
        if(run_lc && loop_start_step != 0 && loop_start_step + lc_period <= step) {
            cout << "lc_running" << endl;
          lc_running = true;
          if(new_loop) {
            loop_start_step = step;
          }
          else {
            loop_start_step = 0;
          }
          new_loop = false;
        }
        else {
          lc_running = false;
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
            isam2.update(newFactors, newVariables, params);

            NonlinearFactorGraph dummy_nfg;
            Values dummy_vals;
            int iter = 0;

            if(run_lc) {
              while(1) {
                isam2.update(dummy_nfg, dummy_vals);
                Values estimate = isam2.calculateEstimate();

                double chi2 = chi2_red(isam2.getFactorsUnsafe(), estimate);

                cout << "iter = " << iter << ", chi2 = " << chi2
                  << ", graph_error = " << isam2.getFactorsUnsafe().error(estimate) << endl;

                if(abs(chi2) < epsilon) {
                  break;
                }

                if(abs(last_chi2 - chi2) < d_error) {
                  break;
                }

                last_chi2 = chi2;

                iter++;

                if(iter >= max_iter) {
                  break;
                } 
              }
            }


            if(lc_running) {
              cout << "update lc" << endl;
              vio_isam2->update(vioNewFactors, vioNewVariables, params);
              Values vio_estimate = vio_isam2->calculateEstimate();
              Values lc_estimate = isam2.calculateEstimate();

              for(int i = loop_start_step; i > 0 && i < step; i++) {
                Pose lcPrevPose = lc_estimate.at<Pose3>(i - 1);
                Pose vioPrevPose = vio_estimate.at<Pose3>(i - 1);
                Pose vioCurPose = vio_estimate.at<Pose3>(i);
                Pose newPose = lcPrevPose * vioPrevPose.inverse() * vioCurPose;
                lc_estimate.update(i, newPose);
              }
              vio_estimate = lc_estimate;
              vio_isam2 = std::make_shared<ISAM2>(isam2params);
              vio_isam2->update(isam2.getFactorsUnsafe(), vio_estimate, params);
            
            }
            else {
              cout << "update vio" << endl;
              vio_isam2->update(vioNewFactors, vioNewVariables, params);
            }
              cout << "end vio update" << endl;

            if(step % print_frequency == 0) {
                // estimate = isam2.calculateEstimate();
                // cout << "Theta = " << endl;
                // estimate.print();
            }
            if(step >= num_steps) {
                break;
            }

            Values vio_estimate = vio_isam2->calculateEstimate();
            Values lc_estimate2 = isam2.calculateEstimate();
            
            double chi2 = chi2_red(isam2.getFactorsUnsafe(), vio_estimate);
            double lc_chi2 = chi2_red(isam2.getFactorsUnsafe(), lc_estimate2);
            cout << "chi2: " << chi2 << " lc_chi2: " << lc_chi2 << endl;

            newVariables.clear();
            newFactors = NonlinearFactorGraph();

            vioNewVariables.clear();
            vioNewFactors = NonlinearFactorGraph();

            string outfile = dataset_outdir + "/step-" + to_string(step) + "_traj.txt";
            ofstream fout(outfile);

            if(!fout.is_open()) {
              cerr << "Cannot open file: " << outfile << endl;
              exit(1);
            }

            vio_estimate.print_kitti_pose3(fout);

            string lc_outfile = dataset_outdir + "/step-" + to_string(step) + "_lc-traj.txt";
            ofstream lc_fout(lc_outfile);

            if(!lc_fout.is_open()) {
              cerr << "Cannot open file: " << lc_outfile << endl;
              exit(1);
            }

            Values lc_estimate = isam2.calculateEstimate();
            lc_estimate.print_kitti_pose3(lc_fout);


        }
        K_count++;

    }



}
