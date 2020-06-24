#ifndef WOODSEER_BRANCH_H
#define WOODSEER_BRANCH_H

#include <stdlib.h>
#include <math.h>

namespace WoodSeer {
    class Branch {
        protected:
            double z_start, z_end;
            double theta_center, sigma_theta_min, sigma_theta_rate;
            double sigma_z_min, sigma_z_rate;

            double surface_amplitude;
            double density_amplitude;
            double alpha1 ;
            double alpha2  ;
            double highFreq ;
            double alpha1c ;
            double alpha2c  ;
            double AmpLowFreqc ;
            double lowFreq ;
            double AmpHighFreq1 ;
            double AmpHighFreq2;
            double AmpLowFreq;
            double e;
            double amplRoot;
            double R;
            double tau;
            double rc;
            double R1;
            double rc_z_min;
            double rc_z_rate;
            double delta;
            double alpha;
        public:
            Branch(bool randomize=true, double theta = 0.0) {
                R=0;
                R1=0.2;
                rc=0;
                surface_amplitude = 1.0;
                density_amplitude = 1.0;
                z_start = 0.0;
                z_end = 0.5;
                theta_center = theta;
                sigma_theta_min = 0.05;
                sigma_theta_rate = 0.05;
                sigma_z_min = 0.05;
                sigma_z_rate = 0.2;
                rc_z_min = 0.1;
                rc_z_rate = 0.2;
                if (randomize) {
                    delta= 2*random()%2-1;
		    alpha1 = -drand48()*0.002+0.4;
		    alpha2 = -drand48()*0.002+0.1;
		    alpha1c = -drand48()*0.002+0.45;
		    alpha2c = -drand48()*0.002+0.15;
	            highFreq = 80+10*drand48();
		    lowFreq = 1+2*drand48();
		    AmpHighFreq1 = drand48()*0.02;
		    AmpHighFreq2 = drand48()*0.02;
		    AmpLowFreq = 0.002+drand48()*0.01;
		    AmpLowFreqc = 0.02+drand48()*0.01;
		    e = 0.2-drand48()*0.04;
		    tau=0.1;
		    amplRoot = 0.1+drand48()*0.001;
            theta_center = theta + randt()/3;
            z_start = drand48()*0.01;
            z_end = z_start + (random()%8)*0.1;
            sigma_theta_rate = 0.02 + drand48()*0.2;
            sigma_z_rate = 0.02 + drand48()*0.4;
            alpha = 0.1*(random()%10);
            
                }
            }

            double getMaxDensity() const {
                return density_amplitude;
            }

            double getMaxSurface(double r) const {
                double sigma_theta = sigma_theta_min+r*sigma_theta_rate;
                double sigma_z = sigma_z_min+r*sigma_z_rate;
                return r * (1 + surface_amplitude * hypot(sigma_theta,sigma_z/r));
            }

            double density(double r, double theta, double z) {
                double texture= AmpHighFreq1*(cos(highFreq*theta)+sin(highFreq*theta))+AmpHighFreq2*(cos(highFreq*z)+sin(highFreq*z));
                double shape=3*AmpLowFreqc*(cos(lowFreq*theta)+sin(lowFreq*theta)) + 3*alpha1c/(1-e*cos(theta))+3*alpha2c/(1-e*cos(z));
                double root_shape=amplRoot*exp(-z/10);
                R = r*(0.5+texture+shape+root_shape);
                double mu_theta = theta_center;
                double sigma_theta = 0.5*(sigma_theta_min+sigma_theta_rate*r);
                double mu_z = z_start + (z_end-z_start)*(1-exp(-r/tau));
                double sigma_z = 0.2*sigma_z_min ;
                double dtheta = remainder(theta-mu_theta,2*M_PI)/(sigma_theta);
                double dz = (z - mu_z)/sigma_z;
                double Gaussian = exp(-0.5*(dtheta*dtheta+dz*dz));

                return density_amplitude*(Gaussian);
            }

            double surface(double r, double theta, double z) {
                double mu_theta = theta_center;

                
                double sigma_theta = sigma_theta_min+r*sigma_theta_rate;
                double mu_z = z_start + (z_end-z_start)*sqrt(r/0.5);
                double sigma_z = sigma_z_min+r*sigma_z_rate;
                double dtheta = remainder(theta-mu_theta,2*M_PI)/(sigma_theta);
                double dz = (z - mu_z)/sigma_z;
                double texture= AmpHighFreq1*(cos(highFreq*theta)+sin(highFreq*theta))+AmpHighFreq2*(cos(highFreq*z)+sin(highFreq*z));
                double shape=2*AmpLowFreqc*(cos(lowFreq*theta)+sin(lowFreq*theta)) + alpha1c/(1-e*cos(theta))+alpha2c/(1-e*cos(z));
                double root_shape=amplRoot*exp(-z/10);
                R = r*(0.2+texture+2*shape+root_shape);
                return R*(1+0.8*surface_amplitude*hypot(sigma_theta*r,sigma_z)
                    *exp(-0.5*(dtheta*dtheta+dz*dz)));
           }
            double circle(double r, double theta, double z) {
                double texture= AmpHighFreq1*(cos(highFreq*theta)+sin(highFreq*theta))+AmpHighFreq2*(cos(highFreq*z)+sin(highFreq*z));
                double shape=2*AmpLowFreqc*(cos(lowFreq*theta)+sin(lowFreq*theta)) + alpha1c/(1-e*cos(theta))+alpha2c/(1-e*cos(z));
                double root_shape=amplRoot*exp(-z/10);
                R = r*(0.2+texture+2*shape+root_shape);
                r=R;


                return R;
            }
            double radius (double r, double theta, double z){
                double mu_theta = theta_center;


                double sigma_theta = 0.7*sigma_theta_min+0.7*r*sigma_theta_rate;
                double mu_z = z_start + (z_end-z_start)*sqrt(r);
                double sigma_z = sigma_z_min+sigma_z_rate*r;
                double dtheta = remainder(theta-mu_theta,2*M_PI)/sigma_theta;
                double dz = (z - mu_z)/(sigma_z);
                double texture= AmpHighFreq1*(cos(highFreq*theta)+sin(highFreq*theta))+AmpHighFreq2*(cos(highFreq*z)+sin(highFreq*z));
                double shape=AmpLowFreqc*(cos(lowFreq*theta)+sin(lowFreq*theta)) + alpha1c/(1-e*cos(theta))+alpha2c/(1-e*cos(z));
                double root_shape=amplRoot*exp(-z/10);
                R = 1+texture+shape+root_shape;

                return R*(1+surface_amplitude*hypot(sigma_theta*r,sigma_z)
                    *exp(-0.5*(dtheta*dtheta+dz*dz)));

            }

            static double randt() {
                double v = 2*drand48()-1;
                if (v<0) {
                    return -v*v;
                } else {
                    return v*v;
                }
            }
    };
};


#endif // WOODSEER_BRANCH_H
