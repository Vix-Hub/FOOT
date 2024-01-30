import pickle

manual_list =  [8, 10, 13, 14, 19, 22, 24, 49, 62, 63, 68, 70, 71, 85, 89, 100, 102, 118, 120, 124, 127, 137, 148, 150, 158, 160, 195, 196, 197, 205, 206, 
212, 213, 221, 226, 227, 235, 238, 240, 250, 253, 262, 272, 296, 294, 299, 301, 305, 314, 316, 317, 318, 321, 325, 330, 333, 334, 337, 348, 349, 
361, 371, 372, 385, 394, 396, 402, 464] #canvas numbers

modifications = [[47, 9448], [272, 58066]]  # [canvasID, mtID1, mtID2, ...] MTs to remove from canvas
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

with open("final_candidates.pkl", "rb") as file:
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


with open("manual_list.pkl", "wb") as file:
    pickle.dump(final_vertices, file)