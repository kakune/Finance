import os
import matplotlib.pyplot as plt
import pandas as pd
import configparser
from py_module.cpp_module import cmake, run_cpp
from py_module.plot_module import graph_plot

dir_path = os.path.dirname(os.path.realpath(__file__))
build_dir = os.path.join(dir_path, "build")
param_file = os.path.join(dir_path, "parameters/parameters.ini")
output_file = os.path.join(dir_path, "build/out/out")
figure_file = os.path.join(dir_path, "out/graph.png")

sections = ["PARAM_BS","PARAM_SABR","PARAM_SV"]

if __name__=="__main__":
    # read params
    params = configparser.ConfigParser()
    params.read(param_file)

    # make csv_file_names from params
    csv_files = [os.path.join(dir_path, "out/result_") + str(params["COMMON"]["name"]) + "_" + str(params[section][params["COMMON"]["name"]]) + ".csv" for section in sections]
    # make labels from params
    labels = [str(params["COMMON"]["name"]) + " = " + str(params[section][params["COMMON"]["name"]])  for section in sections]


    cmake(build_dir)
    
    fig = None
    ax = None
    for section, csv_file, label in zip(sections, csv_files, labels):
        # run C++
        process = run_cpp(output_file = output_file, args = (param_file, section))
        # write the result to a csv file
        with open(csv_file, 'w') as file:
            file.write(process.stdout)
        data_list = pd.read_csv(csv_file, header=0).T.values

        # draw a graph
        fig, ax = graph_plot(data_list[0], data_list[1], fig = fig, ax = ax, label = label)
        # print(process.stdout)

    plt.savefig(figure_file)