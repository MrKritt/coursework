vote_result = [[] for _ in range(len(SVM))]
for i in range(len(ANN)):
    vote_result[i].append(ANN[i][0])
    vote_result[i].append(SVM[i][0])
    vote_result[i].append(KNN[i][0])
    vote_result[i].append(DT[i][0])
    vote_result[i].append(RF[i][0])

for i in range(len(ANN)):
    res = np.mean(vote_result[i])
    if res >= 0.8:
        print('Sequence ', i + 1, ' is a probable Type 6 Effector')
    else:
        print('Sequence ', i + 1, ' is not a Type 6 Effector')