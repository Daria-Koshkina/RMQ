#include <iostream>
#include <vector>
#include <cmath>

const int64_t mod = pow(2, 32);

struct Node {
    int64_t first;
    int64_t second;
};

class RMQ
{
public:
    int block_size, block_count;
    std::vector<int64_t> array;
    std::vector<std::vector<int64_t>> sparse_table;
    std::vector<std::vector<std::vector<int>>> blocks;
    std::vector<int> block_types;

    explicit RMQ(std::vector<int64_t> array) {
        this->array = array;
        precompute_lca();
    }

    int get_min(int i, int j) {
        return array[i] < array[j] ? i : j;
    }

    void precompute_lca() {
        // precompute all log values
        int n = array.size();
        block_size = std::max((int64_t)1, (int64_t)log2(n) / 2);
        block_count = (n + block_size - 1) / block_size;

        // precompute minimum of each block and build sparse table
        sparse_table.assign(block_count,
            std::vector<int64_t>((int64_t)log2(block_count) + 1));
        for (int i = 0, j = 0, b = 0; i < n; i++, j++) {
            if (j == block_size) {
                j = 0, b++;
            }
            if (j == 0 || get_min(i, sparse_table[b][0]) == i) {
                sparse_table[b][0] = i;
            }
        }
        for (int64_t l = 1; l <= (int64_t)log2(block_count); l++) {
            for (int64_t i = 0; i < block_count; i++) {
                int64_t ni = i + ((int64_t)1 << (l - 1));
                if (ni >= block_count) {
                    sparse_table[i][l] = sparse_table[i][l - 1];
                }
                else {
                    sparse_table[i][l] =
                        get_min(sparse_table[i][l - 1], sparse_table[ni][l - 1]);
                }
            }
        }

        // precompute mask for each block
        block_types.assign(block_count, 0);
        for (int i = 0, j = 0, b = 0; i < n; i++, j++) {
            if (j == block_size) {
                j = 0, b++;
            }
            if (j > 0 && (i >= n || get_min(i - 1, i) == i - 1)) {
                block_types[b] += 1 << (j - 1);
            }
        }

        // precompute RMQ for each unique block
        int possibilities = 1 << (block_size - 1);
        blocks.resize(possibilities);
        for (int b = 0; b < block_count; b++) {
            int mask = block_types[b];
            if (!blocks[mask].empty()) {
                continue;
            }
            blocks[mask].assign(block_size, std::vector<int>(block_size));
            for (int l = 0; l < block_size; l++) {
                blocks[mask][l][l] = l;
                for (int r = l + 1; r < block_size; r++) {
                    blocks[mask][l][r] = blocks[mask][l][r - 1];
                    if (b * block_size + r < n) {
                        blocks[mask][l][r] =
                            get_min(b * block_size + blocks[mask][l][r],
                                b * block_size + r) - b * block_size;
                    }
                }
            }
        }
    }

    int rmq_in_block(int b, int l, int r) {
        return blocks[block_types[b]][l][r] + b * block_size;
    }

    int64_t rmq(int l, int r) {
        if (l > r) {
            std::swap(l, r);
        }
        int bl = l / block_size;
        int br = r / block_size;
        if (bl == br) {
            int answ = rmq_in_block(bl, l % block_size, r % block_size);
            return array[answ];
        }
        int ans1 = rmq_in_block(bl, l % block_size, block_size - 1);
        int ans2 = rmq_in_block(br, 0, r % block_size);
        int ans = get_min(ans1, ans2);
        if (bl + 1 < br) {
            int l = (int64_t)log2(br - bl - 1);
            int ans3 = sparse_table[bl + 1][l];
            int ans4 = sparse_table[br - (1 << l)][l];
            ans = get_min(ans, get_min(ans3, ans4));
        }
        return array[ans];
    }
};

int64_t get_sum(int64_t n, int64_t m, int64_t a, int64_t b) {
    std::vector<int64_t> array(0);
    std::vector<Node> queries(0);
    for (int64_t i = 1; i <= n; i++) {
        int64_t temp = (i * a + b) % mod;
        array.push_back(temp);
    }
    for (int64_t i = n + 1; i <= 2 * m + n; i += 2) {
        int64_t temp1 = (i * a + b) % mod;
        int64_t temp2 = ((i + 1) * a + b) % mod;
        queries.push_back(Node{ temp1 % n, temp2 % n });
    }
    int64_t sum = 0;
    RMQ test = RMQ(array);
    for (int64_t i = 0; i < queries.size(); i++) {
        int64_t res = 0;
        res = test.rmq(queries[i].first, queries[i].second);
        sum += res;
    }
    return sum;
}

int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    int64_t n, m, a, b;
    while (true)
    {
        std::cin >> n >> m >> a >> b;
        if (n == 0 && m == 0 && a == 0 && b == 0) {
            break;
        }

        int64_t new_result = get_sum(n, m, a, b);
        std::cout << new_result << "\n";
    }

    return 0;
}
