import multiprocessing
import os
import pickle
import time

completed_processes = multiprocessing.Value('i', 0)

def run_instance(viewID):
    os.system(f'python merge_mts_view.py {viewID}')

    with completed_processes.get_lock():
        completed_processes.value += 1

        # Print progress every 10%
        total_processes = len(viewIDs)
        if completed_processes.value % (total_processes // 10) == 0 or completed_processes.value == total_processes:
            elapsed_time = time.time() - start_time
            print(f"Progress: {completed_processes.value}/{total_processes} ({(completed_processes.value / total_processes) * 100:.1f}%), Elapsed Time: {elapsed_time:.2f} seconds")


if __name__ == '__main__':
    # Specify the input and output file names for each instance

    with open("view_info.pkl", "rb") as file:
        data = pickle.load(file)

    Nviews = len(data["actual_views"])
    viewIDs = []
    for i in range(0, Nviews, 200):
        viewIDs.append(i)
    
    inputs = [(viewID,) for viewID in viewIDs]
    start_time = time.time()

    # Start a separate process for each instance
    with multiprocessing.Pool(processes=len(inputs)) as pool:
        pool.starmap(run_instance, inputs)
