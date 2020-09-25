package com.e.minipaint


import java.lang.Math.sqrt
import java.util.*
import kotlin.math.pow



val k = floatArrayOf(0.0F, 0.0F, 0.0F, 0.0F, 1.0F ,2.0F,3.0F,4.0F, 4.0F, 4.0F ,4.0F )
val c_0 = floatArrayOf(0.0F, 0.0F)
val c_1 = floatArrayOf(1.0F, 0.0F)
val c_2 = floatArrayOf(1.0F, 1.0F)
val c_3 = floatArrayOf(1.0F, 0.0F)
val c_4 = floatArrayOf(0.0F, 0.0F)
val c = mutableListOf<FloatArray>(c_0, c_1, c_2, c_3,c_4)

val u = floatArrayOf(0.0F, 1.0F, 2.0F, 3.0F, 4.0F)


fun main(){
    val S = Spline(k,c, u, 3)
    println(S(2.0F))
}


class BezierInterpolate() {
    /*
    *TO DO:
    * - Make sure solveTridiagonal works ---   Check!
    * - Store the control points so that the curve can be evaluated
    * - Implement another class that uses this class to find the bezier fit with the minimum points
    *      needed to satisfy error less then some tol.
    * - The above would be fitCurve from TimeParametrizedSpline, ported to Kotlin and using this
    *      Bezier implementation instead of B-splines.
    * - Maybe this should not be a class, but rather have a function that returns an object "Bezier",
    *        which has a call method and the control points as property
    * -
     * */

    fun getBezierCoef(x: FloatArray, y: FloatArray): List<FloatArray>{
        /*
        * Finds the handlepoints given a set of points to interpolate.
        *
        * Has to handle x and y separately as long as solveTridiagonal only handles vectors as b
        *
        * Parameters:
        *  x - x-values of the points to interpolate
        *  y - y-values of the points to interpolate
        * Returns:
        *  listof(
        *   A_x: The x_values of handle point A
        *   A_y: The y_values of handle point A
        *   B_x: The x_values of handle point B
        *   B_y: The y_values of handle point B
        *  )
         */
        val n = x.size - 1
        val C = listOf(FloatArray(n-1, { i -> 1F + (i == n-2).float}),
            FloatArray(n, { i -> 4F - 2 * (i == 0).float + 3 * (i == (n-1)).float}),
            FloatArray(n- 1, {i -> 1F})
        )
        val P_x = FloatArray(n, { i -> 2 * (2 * x[i] + x[i + 1])})
        P_x[0] = x[0] + 2 * x[1]
        P_x[n-1] = 8 * x[n - 1] + x[n]
        val P_y = FloatArray(n, { i -> 2 * (2 * y[i] + y[i + 1])})
        P_y[0] = y[0] + 2 * y[1]
        P_y[n-1] = 8 * y[n - 1] + y[n]
        val A_y = solveTridiagonal(C, P_y)
        val A_x = solveTridiagonal(C, P_x)
        val B_x = FloatArray(n, { i -> if (i < n-1) 2 * x[i + 1] - A_x[i + 1] else (A_x[n - 1] + x[n]) / 2})
        val B_y = FloatArray(n, { i -> if (i < n-1) 2 * y[i + 1] - A_y[i + 1] else (A_y[n - 1] + y[n]) / 2})
        return listOf(A_x, A_y, B_x, B_y)
    }
}

class Spline(val k: FloatArray,val c: MutableList<FloatArray>, val u : FloatArray, val p: Int) {
    /* Class for implementing a cubic spline in B-spline form, see https://en.wikipedia.org/wiki/B-spline.
     A cubic spline is a piecewise polynomial, with  continuous second derivatives. The spline S is
     represented as a linear combination of cardinal B-spline basis functions,
     S(u) = sum_i=1^n c_i B_i(u).
     The parameter u is the normalized arc length, u in [0,1]. The knot vector k is not necessarily
     equidistantly spaced, i.e. k_i - k_i+1 may vary with i.
     */

    fun interpolate(u: FloatArray, x: Array<FloatArray>) {
        return
    }

    fun deBoor(u: Float, index: Int, k: FloatArray, c: List<FloatArray>, p: Int): FloatArray {
        /*
        Evaluate S(x). Implementation of deBoors algorithm, see
        https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
        @param x : Double. Point where to evaluate S(x)
        @param index:  The knot interval where x belongs  k[i] <= x < k[i+1]
        @param k : List<Double>. Knot vector. Should be padded in end points
        @param c : List<DoubleArray>. Control point vector
         */

        println(index)
        val d = MutableList(p + 2) { j -> c[j + index - 3 - p] }
        printA(d)
        var alpha: Float = 0.0F
        for (r in 1..p) {
            for (j in p + 1 downTo r) {
                println("r: $r, j: $j, index: $index")
                alpha = (u - k[j + index - p]) / (k[j + 1 + index - r] - k[j + index - p])
                d[j] = (1.0F - alpha) * d[j - 1] + alpha * d[j]
            }
        }

        return d[p]

    }

    operator fun invoke(u: Float): FloatArray {
        val index = bin(u, this.k)
        return deBoor(u, index, this.k, this.c, this.p)

    }


}

