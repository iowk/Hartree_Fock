#include <bits/stdc++.h>
#include "test_classes.h"
#include "../class/class.h"
#include "../src/utils.h"

using namespace std;

void testPrimitive(){
    Primitive p1(3.0, -0.5);
    Primitive p2(2.0, 1.5);
    assert(isEqual(p1.alpha, 3.0) && "Incorrect alpha value for primitive p1");
    assert(isEqual(p1.c, -0.5) && "Incorrect coefficient value for primitive p1");
    assert(isEqual(p2.alpha, 2.0) && "Incorrect alpha value for primitive p2");
    assert(isEqual(p2.c, 1.5) && "Incorrect coefficient value for primitive p2");
    assert(isEqual(getZeta(p1, p2), 5.0) && "Incorrect zeta value between p1 and p2");
    assert(isEqual(getXi(p1, p2), 1.2) && "Incorrect xi value between p1 and p2");
    cout << "Test primitive passed" << endl;
    return;
}
void testBasis(){
    Primitive p1(3.0, -0.5);
    Primitive p2(2.0, 1.5);
    vector<Primitive> primitives;
    primitives.push_back(p1);
    primitives.push_back(p2);
    Basis b1(primitives);
    assert(isEqual(b1.primitives[0].alpha, 3.0) && "Incorrect alpha value for the first primitive");
    assert(isEqual(b1.primitives[0].c, -0.5) && "Incorrect coefficient value for the first primitive");
    assert(isEqual(b1.primitives[1].alpha, 2.0) && "Incorrect alpha value for the second primitive");
    assert(isEqual(b1.primitives[1].c, 1.5) && "Incorrect coefficient value for the second primitive");
    assert(b1.n_primitive == 2 && "Incorrect n_primtive value for basis b1");
    primitives.clear();
    cout << "Test basis passed" << endl;
    return;
}