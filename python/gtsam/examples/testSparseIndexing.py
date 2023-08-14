from scipy.sparse import csc_matrix, csr_matrix

if __name__ == "__main__":
    data = [13, 24, 35, 78, 65]
    rows = [1, 2, 3, 7, 6]
    cols = [3, 4, 5, 8, 5]

    m = csc_matrix((data, (rows, cols)), shape=(10, 10))
    print(m)

    r = [1, 6, 3, 4]
    c = [3, 5, 5, 4]

    indices = [(1, 3), (6, 5), (3, 5), (4, 4)]

    m[r, c] = [1, 1, 1, 1]

    print(m)
    print(m.nonzero())

