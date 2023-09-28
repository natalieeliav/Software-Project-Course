import sys
import sklearn.metrics
from sklearn.metrics import silhouette_score
import symnmf

epsilon = 0.001
all_vectors = []
centers = []
d = {}
default_iter = 200
muK = []


def init(k):
    for i in range(k):
        centers.append(all_vectors[i])
        d[i] = []
        muK.append(10)

def import_points(data):
    for line in data:
        vector = line.split(",")
        for i in range(len(vector)):
            vector[i] = float(vector[i])
        all_vectors.append(vector)

def distance(vector1, vector2):
    sum = 0
    for i in range(len(vector1)):
        sum += (vector1[i]-vector2[i])**2
    return sum**0.5

def choose_center(vector):
    min = distance(vector, centers[0])
    index_chosen = 0
    for i in range(1, len(centers)):
        d = distance(vector, centers[i])
        if d<min:
            min = d
            index_chosen = i
    return index_chosen


def update_center(index):
    last_center = centers[index]
    vectors = d[index]
    vector_size = len(vectors[0])
    new_center = [0 for i in range(vector_size)]
    num_of_vectors = len(vectors)
    for i in range(vector_size):## Lekol i in len(vector)
        for vec in vectors:
            new_center[i] += vec[i]
        new_center[i] = new_center[i]/num_of_vectors
    centers[index] = new_center
    muK[index] = distance(last_center, new_center)

def update_all_centers(k):
    for i in range(k):
        update_center(i)

def is_valid(k, iter):
    if (iter >= 1000) or (iter <= 1) or (type(iter) != int):
        print("Invalid maximum iteration!")
        return False
    if (k <= 1) or (k >= len(all_vectors)) or (type(k) != int):
        print("Invalid number of clusters!")
        return False
    return True    

def is_valid_int(str):
    for char in str:
        if ord(char) < ord('0') or ord(char) > ord('9'):
            return False
    return True

def kMean(K,iter=default_iter):
    for i in range(iter):
        counter = 0
        for i in range(K):
            d[i] = []
        for vector in all_vectors:
            d[choose_center(vector)].append(vector)
        update_all_centers(K)
        for mu in muK:
            if mu < epsilon:
                counter += 1
            else :
                break
        if counter == K:
            break

def get_Label(clusters) :
    label = []
    for vector in all_vectors :
        label.append(choose_center(vector))
    return label

def H_centers(H) :
    return [row.index(max(row)) for row in H] 

def main():
    args_amount = len(sys.argv)
    if args_amount != 3:
        print("An Error Has Occurred")
        return
    else:
        file_in = sys.argv[2]
        file = open(file_in)
        import_points(file)
        str_k = sys.argv[1]
        if is_valid_int(str_k):
            k = int(str_k)
        else:
            print("Invalid number of clusters!")
            return
        if not is_valid(k, default_iter):
            return
        init(k)
        kMean(k)
        h = symnmf.final_H(all_vectors, len(all_vectors), len(all_vectors[0]), k)
        kmeans_labels = get_Label(centers)
        symnmf_labels = H_centers(h)
        silhouette_kmeans = silhouette_score(all_vectors, kmeans_labels)
        silhouette_symnmf = silhouette_score(all_vectors, symnmf_labels)
        print("nmf:", silhouette_symnmf)
        print("kmeans:", silhouette_kmeans)
        return
    

if __name__ == "__main__":
    main()





