#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, mdivideLeftSpd) {
  auto f = [](const auto& x, const auto& y) {
    if (x.rows() != x.cols())
      return stan::math::mdivide_left_spd(x, y);
    auto x_sym = ((x + x.transpose()) * 0.5).eval();  // sym for finite diffs
    return stan::math::mdivide_left_spd(x_sym, y);
  };

  // signature 1 of 2: matrix-matrix
  Eigen::MatrixXd aa(1, 1);
  aa << 1;
  Eigen::MatrixXd bb(1, 1);
  bb << 2;
  stan::test::expect_ad(f, aa, bb);

  // signature 2 of 2: matrix-vector
  Eigen::VectorXd cc(1);
  cc << 3;
  stan::test::expect_ad(f, aa, cc);

  Eigen::MatrixXd m00(0, 0);
  Eigen::VectorXd v0(0);
  stan::test::expect_ad(f, m00, v0);
  stan::test::expect_ad(f, m00, m00);

  Eigen::MatrixXd a(2, 2);
  a << 2, 3, 3, 7;

  Eigen::MatrixXd b(2, 2);
  b << 3, 0, 0, 4;

  Eigen::MatrixXd c(2, 2);
  c << 12, 13, 15, 17;

  Eigen::VectorXd d(2);
  d << 12, 13;

  // matrix, matrix : these have spd first args as required
  stan::test::expect_ad(f, a, a);
  stan::test::expect_ad(f, b, b);
  stan::test::expect_ad(f, a, c);
  stan::test::expect_ad(f, b, c);
  // matrix, vector : ditto
  stan::test::expect_ad(f, a, d);
  stan::test::expect_ad(f, b, d);
  stan::test::expect_ad(f, a, d);
  stan::test::expect_ad(f, b, d);

  Eigen::MatrixXd m33 = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd m44 = Eigen::MatrixXd::Zero(4, 4);
  Eigen::VectorXd v3 = Eigen::VectorXd::Zero(3);
  Eigen::VectorXd v4 = Eigen::VectorXd::Zero(4);
  Eigen::RowVectorXd rv3 = Eigen::RowVectorXd::Zero(3);
  Eigen::RowVectorXd rv4 = Eigen::RowVectorXd::Zero(4);

  // exceptions: not symmetric
  stan::test::expect_ad(f, c, a);
  stan::test::expect_ad(f, c, d);

  // FIXME(carpenter): double throws (correct) but var doesn't (incorrect)
  // exceptions: not pos def
  // stan::test::expect_ad(f, m33, m33);
  // stan::test::expect_ad(f, m33, v3);

  // exceptions: wrong sizes
  stan::test::expect_ad(f, m33, m44);
  stan::test::expect_ad(f, m33, v4);

  // exceptions: wrong types
  stan::test::expect_ad(f, m33, rv3);
}
