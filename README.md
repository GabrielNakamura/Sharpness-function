# Sharpness-function
Function to test the significance of partitions in cluster analysis of sampling units. The method is described in [Pillar, 1999](https://esajournals.onlinelibrary.wiley.com/doi/pdf/10.1890/0012-9658%281999%29080%5B2508%3AHSAC%5D2.0.CO%3B2).

#arguments
inputs:
x= a matrix object with sampling units in rows and descriptors in columns;

k= a scalar. Value greater than 1 indicating the number of groups to be tested regarding the significance;

method= character, indicating which measure will be applied in x to calculate the the distance matrix used in cluster analysis. Can be one of the methods presented in [vegdist](https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/vegdist) function, default is "euclidean";

iterations= scalar, indicating the number of bootstrap iterations used in the procedure to calculate the significance. Default is 999.

output:
matrix with five columns containing the values of mean observed statistic, mean null statistic, their respectives standard deviation and the p value associated with the observed statistic.
