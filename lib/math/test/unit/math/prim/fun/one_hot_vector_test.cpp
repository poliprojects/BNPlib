#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/expect_matrix_eq.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, one_hot_vector) {
  using Eigen::VectorXd;
  using stan::math::one_hot_vector;

  for (int K = 1; K < 5; K++) {
    for (int k = 1; k <= K; k++) {
      Eigen::VectorXd y = Eigen::VectorXd::Zero(K);
      y[k - 1] = 1;
      expect_matrix_eq(y, one_hot_vector(K, k));
    }
  }
}

TEST(MathFunctions, one_hot_vector_throw) {
  using stan::math::one_hot_vector;
  int K = 5;
  int k = 2;

  EXPECT_THROW(one_hot_vector(-1, k), std::domain_error);
  EXPECT_THROW(one_hot_vector(0, k), std::domain_error);
  EXPECT_THROW(one_hot_vector(K, K + 1), std::domain_error);
  EXPECT_THROW(one_hot_vector(K, 0), std::domain_error);
  EXPECT_THROW(one_hot_vector(K, -1), std::domain_error);
}
