import numpy as np
from scipy.special import softmax


input = np.array([[1,0],[0,1],[1,1],[0,0]])

projK = np.array([[1, 0], [0, 0]])
projQ = np.array([[0, 0], [1, 0]])
projV = np.array([[1, 0], [0, 1]])
biasK = np.array([0, 0])
biasQ = np.array([0, 0])
biasV = np.array([0, 0])

Keys1 = np.dot(input, projK[0]) + biasK[0]

Keys2 = np.dot(input, projK[1]) + biasK[1]

Queries1 = np.dot(input, projQ[0]) + biasQ[0]

Queries2 = np.dot(input, projQ[1]) + biasQ[1]

Values1 = np.dot(input, projV[0]) + biasV

Values2 = np.dot(input, projV[1]) + biasV

Logits1 = np.dot(Queries1, Keys1.T)

Logits2 = np.dot(Queries2, Keys2.T)

AttScores1 = softmax(Logits1, axis=1)

AttScores2 = softmax(Logits2, axis=1)

Result1 = np.dot(AttScores1,Values1[:,0])

Result2 = np.dot(AttScores2, Values2[:,0])

result = np.stack((Result1,Result2)).T

print(result)