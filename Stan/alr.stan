  #include util.stan

  /* additive log ratio transformation
   *
   * Args:
   *   x: a simplex (but cannot declare argument as type simplex)
   *   i: an integer giving the index of the reference element
   * Returns:
   *   the transformed vector with 1 element fewer than x
   */
  vector alr(vector x, int i) {
    int l = num_elements(x);

    if (l < i) {
      reject("x must be at least as long as i");
    }

    if (l < 2) {
      reject("x must have at least 2 elements");
    }

    if (abs(sum(x) - 1) > 1e-6) {
      reject("x must be a simplex, found ", sum(x));
    }

    if (x[i] <= 0) {
      reject("the reference element of x must be positive, found ", x[i]);
    }

    array[l - 1] int idx; 

    if (i == 1) {
      idx = sequence(2, l);
    } else if (i == l) {
      idx = sequence(1, l - 1);
    } else {
      idx = append_array(sequence(1, i - 1), sequence(i + 1, l));
    }

    return log(x[idx]) - log(x[i]);
  }

  /* inverse additive log ratio transformation
   *
   * Args:
   *   x: a vector
   *   i: an integer giving the index that was used for the ALR transform
   * Returns:
   *   the inverse of alr(x)
   */
  vector alrinv(vector x, int i) {
    int l = num_elements(x);

    if (l < i - 1) {
      reject("x must be at least as long as i - 1");
    }

    vector[l + 1] res;

    if (i == 1) {
      res = append_row(0, x);
    } else if (i == l + 1) {
      res = append_row(x, 0);
    } else {
      res = append_row(append_row(x[:(i - 1)], 0), x[i:]);
    }

    return softmax(res);
  }
