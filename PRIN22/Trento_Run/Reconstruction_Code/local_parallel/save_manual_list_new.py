import pickle

manual_list =  [] #canvas numbers

manual_list_cuts = [ 1, 4, 5, 6, 8, 13, 14, 15, 17, 21, 25, 26, 27, 28, 29, 33, 34, 35, 36, 43, 45, 47, 55, 57, 58, 59, 61, 62, 
66, 67, 68, 71, 72, 78, 80, 81, 82, 87, 88, 90, 92, 95, 98, 101, 102, 103, 106, 107, 114, 116, 117, 118, 119, 
123, 124, 125, 126, 130, 134, 136, 137, 139, 140, 142, 147, 155, 158] #cuts

CUTS = 1
filename = "final_candidates.pkl"
outname = "manual_list.pkl"

if (CUTS):
    manual_list = manual_list_cuts
    filename = "final_candidates_cuts.pkl"
    outname = "manual_list_cuts.pkl"

modifications = []  # [canvasID, mtID1, mtID2, ...] MTs to remove from canvas
to_be_modded, to_remove2 = [], [] 

for entry in modifications:
    to_be_modded.append(entry[0])
    to_remove2.append([])
    for i in range(1, len(entry)):
        to_remove2[-1].append(entry[i])


to_remove = []
print(len(manual_list))
counter = 0
for entry in to_remove:
    manual_list.remove(manual_list[entry-counter])
    counter += 1

with open(filename, "rb") as file:
    final_candidates = pickle.load(file)

final_vertices = []
duplicates_check = []
for canvasID in manual_list:
    if (not (canvasID in to_be_modded)):
        for entry in final_candidates[canvasID]:
            if (not (entry in duplicates_check)):
                duplicates_check.append(entry)
            else:
                continue
        final_vertices.append(final_candidates[canvasID])
    else:
        temp = final_candidates[canvasID]
        for to_be_removed in to_remove2[to_be_modded.index(canvasID)]:
            print(temp, to_be_removed, final_candidates[canvasID-1], final_candidates[canvasID+1])
            temp.remove(to_be_removed)
        if (not (entry in duplicates_check)):
                duplicates_check.append(entry)
        else:
            continue
        final_vertices.append(final_candidates[canvasID])


with open(outname, "wb") as file:
    pickle.dump(final_vertices, file)
