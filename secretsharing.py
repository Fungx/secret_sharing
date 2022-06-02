from galois import GF


class Shamir:
    """
    Shamir implements Shamir's secret sharing scheme.

    Reference:
    Shamir, Adi. "How to share a secret." Communications of the ACM 22.11 (1979): 612-613.
    """

    def __init__(self, order):
        """
        Initializes with the assigned finite field.
        :param order: order of the finite field. The order must be a prime power.
        """
        self.__gf = GF(order)

    def share(self, k: int, n: int, secret: int) -> list:
        """
        Splits a secret into `n` shares using the (k,n)-threshold scheme.
        Any subset of `k` of these shares can reveal the secret.
        :param secret: the shared secret.
        :param k: threshold of the scheme.
        :param n: number of shares.
        :return: list of share tuples, e.g., `[(1, 185), (2, 101), (3, 218)]`.
                 Supposes the generated polynomial is `p(x)`,
                 the first element of each tuple is `x` and the second is the value of `p(x)`.
        """
        coeffs = self.__gf.Random(k)
        coeffs[k - 1] = self.__gf(secret)
        return [(x, int(self.__poly(x, coeffs))) for x in range(1, n + 1)]

    def reveal(self, shares: list) -> int:
        """
        Reveals the shared secret from the given shares using the Lagrange interpolation.
        A wrong secret may be returned if the size of given shares is less than `k`.
        :param shares: list of share tuples, e.g., `[(1, 185), (2, 101), (3, 218)]`.
        :return: the shared secret.
        """
        x, y = zip(*self.__gf(shares))
        product = self.__gf(1)
        # precomputes the numerator in li
        for xi in x:
            product *= xi

        sums = self.__gf(0)
        for i in range(len(x)):
            denominator = self.__gf(1)
            for j in range(len(x)):
                if i == j:
                    continue
                denominator *= x[i] - x[j]
            li = (product / x[i]) / denominator  # li(0)
            sums += y[i] * li
        return int(sums)

    def __poly(self, x, coeffs):
        """
        Generates a `y` on the graph of given polynomial in `x`.
        The polynomial is given by the list of `coeffs`.
        """
        y = self.__gf(0)
        for index, val in enumerate(coeffs[::-1]):
            y += x ** index * val
        return y


if __name__ == '__main__':
    # (3,5) sharing scheme
    k, n = 3, 5
    secret = 150
    print(f'Original Secret: {secret}')
    ss = Shamir(256)
    # Generation of shares
    shares1 = ss.share(k, n, secret)
    print(f'Shares1: {shares1}')
    shares2 = ss.share(k, n, secret)
    print(f'Shares2: {shares2}')

    # Secret Reconstruction
    pool = shares1[0:3]
    print(f'Reconstruct with Shares1 {pool}. Reconstructed Secret: {ss.reveal(pool)}')
    pool = shares2[0:3]
    print(f'Reconstruct with Shares2 {pool}. Reconstructed Secret: {ss.reveal(pool)}')
    pool = [(shares1[i][0], ((shares1[i][1] + shares2[i][1]) % 256)) for i in range(3)]
    print(f'Reconstruct with Shares1 + Shares2 {pool}. Reconstructed Secret: {ss.reveal(pool)}')
    pool = [(shares1[i][0], ((shares1[i][1] + 80) % 256)) for i in range(3)]
    print(f'Reconstruct with Shares1 + 80 {pool}. Reconstructed Secret: {ss.reveal(pool)}')
    pool = [(shares1[i][0], ((shares1[i][1] * 3) % 256)) for i in range(3)]
    print(f'Reconstruct with Shares1 * 3 {pool}. Reconstructed Secret: {ss.reveal(pool)}')
