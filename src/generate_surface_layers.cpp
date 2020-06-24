#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>


#include "Log.h"
#include "Branch.h"
using namespace WoodSeer;


int main(int argc, char *argv[]) {
    int num = 1;
    if (argc > 1) {
        num = atoi(argv[1]);
    }

    unsigned int N = 256;
    unsigned int N_z = 32;
    double scale = 1.4/N;
    double scale_z = 0.7/N_z;
    srand48(time(NULL));

    double s = 0.0;
    double s_comp = 0.0;
    double upper_thres = 0.01;
    double lower_thres = -3e-2;
    double upper_thres_cs = 0.02;
    double lower_thres_cs = -0.04;
    double upper_thres_complete = 0.0;
    double lower_thres_complete = -2e-2;
    while (num > 0) {
        unsigned int nb_branches = 2+random()%2;
        printf("Log num_branches is %03d num is %03d\n",num,nb_branches); 
        Log L(nb_branches);
        for (unsigned int k=0;k<N_z;k++) {
            double z = k*scale_z;
            cv::Mat_<float> surface(N,N,0.0);
            cv::Mat_<float> surface_comp(N,N,0.0);
            cv::Mat_<float> density(N,N,0.0);
            for (unsigned int i=0;i<N;i++) {
                double x = (double(i) - N/2.)*scale;
                for (unsigned int j=0;j<N;j++) {
                    double y = (double(j)-N/2.)*scale;
                    double r = hypot(x,y);
                    double theta = atan2(y,x);
                    double r_s = L.surface(theta,z);
                    double d = L.density(r,theta,z);
                    double c_s = L.circle(theta,z);
                    double output=2*r-r_s;
                    double output_cs=2*r-c_s;
                    if (output<upper_thres && output>lower_thres && not(output_cs<upper_thres_cs && output_cs>lower_thres_cs)){s = 1;}
                    else{s = 0;}
                    if (output<upper_thres_complete && output>lower_thres_complete){s_comp = 1;}
                    else{s_comp = 0;}
                    surface(i,j)=s;
                    surface_comp(i,j)=s_comp;
                    density(i,j)=d;
                    }
                }
            
            cv::Mat_<uint8_t> surf_comp_out(surface_comp.size());
            surface_comp = 255.0*(surface_comp);
            surface_comp.convertTo(surf_comp_out,CV_8U);
            char tmp[128];
            sprintf(tmp,"surface_%03d_comp_%03d.png",num,k);
            cv::imwrite(tmp,surf_comp_out);


            cv::Mat_<uint8_t> surf_out(surface.size()); 
            surface = 255*(surface);
            surface.convertTo(surf_out,CV_8U);
            char surf[128];
            sprintf(surf,"surface_%03d_%03d.png",num,k);
            cv::imwrite(surf,surf_out);
    
    
            cv::Mat_<uint8_t> d_out(density.size()); 
            density = 255*(density);
            density.convertTo(d_out,CV_8U);
            char d[128];
            sprintf(d,"density_%03d_%03d.png",num,k);
            cv::imwrite(d,d_out);
    
    
    
            }


        num -= 1;
    }
    return 0;
};





