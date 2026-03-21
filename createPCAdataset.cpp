#include "PCA.h"

int main() {
    std::string folder = "../data/CubicLatticeKnots/L100/";
    int reductionFactor = 1;
    auto dataset = createDataset(folder, reductionFactor);
    PCA pca;
    pca.fit(dataset);
}