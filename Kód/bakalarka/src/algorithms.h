#ifndef ALGORITHMS_H
#define ALGORITHMS_H

numeric Newton_Raphson(realsymbol my_x, const ex &f, const ex &df,
                       numeric starting_point) {

  numeric precision = 1e-6;

  numeric iter = starting_point;
  int iterations = 0;
  while (abs(f.subs(my_x == iter).evalf()) > precision && iterations < 1'000) {
    iterations++;
    assertm(abs(df.subs(my_x == iter).evalf()) > 10e-6,
            "Division by 0 in N-R method!");
    iter -= ex_to<numeric>(f.subs(my_x == iter).evalf() /
                           df.subs(my_x == iter).evalf());
  }
  // cout << "Number of iterations: " << iterations << endl;

  return iter;
}

#endif
