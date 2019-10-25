import numpy as np
import matplotlib.pyplot as plt


def compress_image(origin_image, rate=0.8, compress_alpha=False):
    result = np.zeros(origin_image.shape)
    for chan in range(4):
        U, sigma, V = np.linalg.svd(origin_image[:, :, chan])
        n_sigmas = 0
        temp = 0
        while (temp / np.sum(sigma)) < rate:
            temp += sigma[n_sigmas]
            n_sigmas += 1
        S = np.zeros((n_sigmas, n_sigmas))
        for i in range(n_sigmas):
            S[i, i] = sigma[i]
        result[:, :, chan] = (U[:, 0:n_sigmas].dot(S)).dot(V[0:n_sigmas, :])
    for i in range(4 if compress_alpha else 3):
        MAX = np.max(result[:, :, i])
        MIN = np.min(result[:, :, i])
        result[:, :, i] = (result[:, :, i] - MIN) / (MAX - MIN)
    result = np.round(result * 255).astype('int')
    plt.subplot(1, 3, 3)
    plt.imshow(result)
    plt.title('compress by python')


def out_img_data(origin_image):
    with open('img_data.txt', 'w') as f:
        f.write("%d %d\n" % (origin_image.shape[0], origin_image.shape[1]))
        for i in range(origin_image.shape[0]):
            for j in range(origin_image.shape[1]):
                for chan in range(4):
                    f.write("%d%s" % (
                        origin_image[i, j, chan], '\n' if j == origin_image.shape[1] - 1 and chan == 3 else ' '))


def read_and_show():
    with open("img.txt", 'r') as f:
        content = [int(i) for i in f.read().split()]
        n, m = content[:2]
    content = np.array(content[2:]).reshape((n, m, 4))
    plt.subplot(1, 3, 2)
    plt.imshow(content)
    plt.title('compress by C')


def main():
    image_path = 'test.jpg'
    origin_image = plt.imread(image_path)
    plt.figure(figsize=(12, 4))
    plt.subplot(1, 3, 1)
    plt.imshow(origin_image)
    plt.title('origin')
    read_and_show()
    compress_image(origin_image, rate=0.8)
    plt.savefig('result.png')
    plt.show()


if __name__ == "__main__":
    main()
    '''
    arr = np.array([
        [5, 3, 2],
        [3, 5, 1],
        [2, 1, 1]
    ])
    U, sigma, V = np.linalg.svd(arr)
    print(U)
    print(sigma)
    print(V)
    '''
