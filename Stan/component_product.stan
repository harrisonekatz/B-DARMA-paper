/* Component-wise product
 * 
 * Args:
 *   x: an array of vectors
 *   y: a vector to multiply against the segments of x
 *   k: an array giving the number of elements of x by component
 * Returns:
 *   an array of vectors as long as x, each element of which contains the
 *   component-wise product of x with y
 */
array[] vector component_product(array[] vector x, vector y, array[] int k) {
  int N = size(x);
  int C = size(k);

  array[N] vector[C] result;

  for (n in 1:N) {
    int i = 1;

    for (c in 1:C) {
      result[n][c] = dot_product(
        segment(x[n], i, k[c]), 
        segment(y, i, k[c])
      );
      i += k[c];
    }
  }

  return result;
}
