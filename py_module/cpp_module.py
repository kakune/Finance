import subprocess
import os

def compile_cpp(
        source_files,
        output_file = "out",
        args = (),
    ):
    # compile C++
    try:
        compile_command = ["icpx", *source_files, "-o", output_file, *args]
        subprocess.run(compile_command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"returncode:{e.returncode}")
        print(e.stderr)
        exit()

def cmake(build_dir):
    os.makedirs(build_dir, exist_ok=True)
    # CMakeを実行してプロジェクトを構成します。
    subprocess.check_call(["cmake", "-S", ".", "-B", build_dir])
    # CMakeを用いてビルドを実行します。
    subprocess.check_call(["cmake", "--build", build_dir])

def run_cpp(
        output_file = "out",
        args = (),
    ):
    # run C++ program
    try:
        process = subprocess.run([output_file, *args], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"returncode:{e.returncode}")
        print(e.stderr)
        exit()
    return process
