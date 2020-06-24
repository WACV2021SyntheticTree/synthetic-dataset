# Synthetic Dataset Generator

- This project is  for surface and internal density log generation. The number of logs and the resolution of the cross-sections have to be changed inside ./src/generate_surface_layers.cpp
- The build_box.py is a python script for external defect automatic bouding box annotation. It generates .xml file for each surface cross section.  
 

### Tech

Synthetic Dataset Generator uses a number of open source projects to work properly:

* [C++] - generated the synthetic dataset!
* [Python] - surface cross-section labeling 

### Installation

Dillinger requires some C++ and Python packages 

Install the the required packages through the following command line.

```sh
$ sudo apt-get install cmake
$ sudo apt-get install libopencv-dev
$ pip install opencv-python
$ pip install numba
$ pip install pandas 
$ pip install tqdm
```
### Development

To generate the synthetic images first Tab:
```sh
$ mkdir build
$ cd ./build
$ cmake ..
$ make 
$ cd ..
$ ./generate_surface_layers
```

To generate the .xml files containing the external defects bounding boxes parameters
```sh
$ mv build_box.py bin
$ cd bin
$ python build_box.py
```
License
----
The MIT License (MIT)


