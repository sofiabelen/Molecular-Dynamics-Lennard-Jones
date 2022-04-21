#include <vector>
#include <stdio.h>
#include <iostream>
#include <random>

#include "vector.hpp"

// 3-dim constructor, all alements are set to 0
Vector::Vector() {

    this->dim = 3;

    coord.resize(this->dim, 0);
}

// n-dim constructor, all alements are set to 0
Vector::Vector(const int &dim) {

    this->dim = dim;

    coord.resize(dim, 0);
}

// n-dim constructor, all alements are set to val
Vector::Vector(const int &dim, const double &val) {

    this->dim = dim;

    coord.resize(dim, val);
}

// 3-dim constructor
Vector::Vector(const int &x, const int &y, const int& z) {

    dim = 3;

    coord.resize(dim);

    coord[0] = x;
    coord[1] = y;
    coord[2] = z;
}

Vector Vector::operator +(const Vector &b) {

    Vector sum(this->dim);

    for (int i = 0; i < this->dim; i++) {
        sum.coord[i] = this->coord[i] + b.coord[i];
    }
    return sum;
}

Vector Vector::operator -(const Vector &b) {

    Vector subtr(this->dim);

    for (int i = 0; i < this->dim; i++) {
        subtr.coord[i] = this->coord[i] - b.coord[i];
    }
    return subtr;
}

Vector Vector::operator *(const Vector &b) {

    Vector dotted(this->dim);

    for (int i = 0; i < this->dim; i++) {
        dotted.coord[i] = this->coord[i] * b.coord[i];
    }
    return dotted;
}

Vector Vector::operator *(const double &scalar) {

    Vector product(this->dim);

    for (int i = 0; i < this->dim; i++) {
        product.coord[i] = this->coord[i] * scalar;
    }
    return product;
}

Vector Vector::operator /(const double &scalar) {

    Vector div(this->dim);

    for (int i = 0; i < this->dim; i++) {
        div.coord[i] = this->coord[i] / scalar;
    }
    return div;
}

Vector& Vector::operator *=(const double &scalar) {

    for (int i = 0; i < this->dim; i++) {
        this->coord[i] *= scalar;
    }
    return *this;
}

Vector& Vector::operator /=(const double &scalar) {

    for (int i = 0; i < this->dim; i++) {
        this->coord[i] /= scalar;
    }
    return *this;
}

Vector& Vector::operator +=(const Vector &b) {
    for (int i = 0; i < dim; i++) {
        coord[i] += b.coord[i];
    }
    return *this;
}

Vector& Vector::operator -=(const Vector &b) {
    for (int i = 0; i < dim; i++) {
        coord[i] -= b.coord[i];
    }
    return *this;
}

// Returns length squared
double Vector::len_sq() {
    double sum = 0;
    for (int i = 0; i < dim; i++) {
        sum += coord[i] * coord[i];
    }
    return sum;
}

// Apply periodic boundary conditions
void Vector::wrap(const double &size) {

    for (int i = 0; i < dim; i++) {

        if (coord[i] < 0) {
            coord[i] += size;
        } else if (coord[i] > size) {
            coord[i] -= size;
        }
    }
}

// Print coordinates to standar output
void Vector::show() {
    for (int i = 0; i < dim; i++) {
        printf("%.2f ", this->coord[i]);
    }
    printf("\n");
}
