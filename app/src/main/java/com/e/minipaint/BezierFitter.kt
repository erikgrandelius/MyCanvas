package com.e.minipaint


import java.lang.Math.sqrt
import java.util.*
import kotlin.math.cos
import kotlin.math.pow

<<<<<<< HEAD


val k = floatArrayOf(0.0F, 0.0F, 0.0F, 0.0F, 1.0F ,2.0F,3.0F,4.0F, 4.0F, 4.0F ,4.0F )
val c_0 = floatArrayOf(0.0F, 0.0F)
val c_1 = floatArrayOf(1.0F, 0.0F)
val c_2 = floatArrayOf(1.0F, 1.0F)
val c_3 = floatArrayOf(1.0F, 0.0F)
val c_4 = floatArrayOf(0.0F, 0.0F)
val c = mutableListOf<FloatArray>(c_0, c_1, c_2, c_3,c_4)

val u = floatArrayOf(0.0F, 1.0F, 2.0F, 3.0F, 4.0F)
=======
/*
Here we collect numpy-like methods for FloatArray.
 */
//In order to do math with the truth-value of an expression. Might be better to use if-statement in the lambda
val Boolean.float get() = if (this) 1F else 0F
fun DoubleArray.toString() = Arrays.toString(this)
fun IntArray.toString() =  Arrays.toString(this)
operator fun FloatArray.times(a: FloatArray) = FloatArray(this.size){ i -> this[i] * a[i]}
operator fun FloatArray.times(a: Float) = FloatArray(this.size){ i -> this[i] * a }
operator fun Float.times(a: FloatArray) = FloatArray(a.size){ i -> a[i] * this }
fun FloatArray.pow(p:Float) = FloatArray(this.size){i -> this[i].pow(p)}
fun FloatArray.pow(n:Int) = FloatArray(this.size){i -> this[i].pow(n)}

fun norm(a: FloatArray, p: Float = 2.0f ): Float{
    return a.pow(p).sum().pow(1/p)
}

fun linspace(a: Float, b: Float, n: Int): FloatArray{
    return FloatArray(n){i -> a * (n - 1 - i ) / (n - 1) + b * i / (n - 1)}
}












fun solveTridiagonal(A: List<FloatArray>, b: FloatArray): FloatArray {
    /*
    * Solve Ax = b, for tridiagonal Matrix A.
    * Implements the Thomas algorithm (more or less Gaussian elimination)
    * for solving the tri-diagonal linear system.
    * https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm.
    * Only stable in the diagonally dominant or positive definite case.
    * Input vectors will be altered.
    *
    *
    * @Param: A : List<DoubleArray>. List of Lower diagonal, diagonal, and upper diagonal:
    *                  [l, d, u]. l and u have size n-1, and d has size n.
              b : DoubleArray.  Right hand vector to solve for. Size is n.
    *
    * Returns: x: DoubleArray
               The solution. Size is n.
    */

    val d = b
    val a = A[0].copyOf()
    val b = A[1].copyOf()
    val c = A[2].copyOf()
    var w: Float
    val n = b.size
    for (i in 1..n - 1) {
        w = a[i - 1] / b[i - 1]
        b[i] = b[i] - w * c[i - 1]
        d[i] = d[i] - w * d[i - 1]
    }
    val x = FloatArray(n)
    x[n - 1] = d[n - 1] / b[n - 1]
    for (i in n - 2 downTo 0) {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]
    }
    return x
}




>>>>>>> e36215d7f928cdf5d229fb91068c00a8c73570e8

val u = linspace(0f, 1f, 3)
val x = linspace(0f, 1f, 3)
var y = FloatArray(3) {i -> cos(x[i])}

<<<<<<< HEAD
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
=======
val b = arrayOf(floatArrayOf( 0f, 0f), floatArrayOf(1f, 1f), floatArrayOf(2f,1f), floatArrayOf(3f,0f))



class Spline( x: FloatArray, y: FloatArray, val k: FloatArray) {

    /* Class for implementing a cubic spline in B-spline form, see https://en.wikipedia.org/wiki/B-spline.
     A cubic spline is a piecewise polynomial, with  continuous second derivatives. The spline S is
     represented as a linear combination of cardinal B-spline basis functions,
     S(u) = sum_i=1^n c_i B_i(u).
     The parameter u is the normalized arc length, u in [0,1]. The knot vector k is not necessarily
     equidistantly spaced, i.e. k_i - k_i+1 may vary with i.
     */
    //

    val b = getBezierCoef(x,y)
>>>>>>> e36215d7f928cdf5d229fb91068c00a8c73570e8

    fun getBezierCoef(x: FloatArray, y: FloatArray): Array<Array<FloatArray>>{
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
        val d = Array<Array<FloatArray>> (n){
            i ->
            arrayOf(floatArrayOf(x[i], y[i]),
            floatArrayOf(A_x[i], A_y[i]),
            floatArrayOf(B_x[i], B_y[i]),
            floatArrayOf(x[i + 1], y[i + 1]))
        }
        return d
    }

<<<<<<< HEAD
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
=======
    fun evalBezier(t: Float , b: Array<FloatArray>): FloatArray{
        /*
         Evaluate cubic Bezier curve t -> B(t) defined by control points b0, b1, b2, b3   and  t in [0,1].
         Implements De_Casteljau's_algorithm, see https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm.
         This is more numerically stable than direct evaluation from coefficients.

         @param b:  Array of coefficients, should have size 4.
         @param t:  Point for evaluate, should belong to the interval [0,1].

         @returns d[0]: Value of B(t).
         */
        var d: Array<FloatArray> = b.copyOf()
        for (i in d.size - 2  downTo 0) {
            for (j in 0 .. i )
                d[j] = (1 - t) * d[j] + t * d[j + 1]
        }
        return d[0]
    }
    fun findKnotIndex(u : Float, k : FloatArray): Int {
>>>>>>> e36215d7f928cdf5d229fb91068c00a8c73570e8
        /*
        For values x and bin k, return Array of indices j such that
        k[J[j]]<= x[J[j]] < k[J[j] + 1 ]
        @param x, array of values
        @param k, ordered array of bins

        returns: J, array of indices
         */

<<<<<<< HEAD
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
=======
        val k_left = k.dropLast(1)
        val k_right = k.drop(1)
        return k_left.indices.filter { (k_left[it] <= u)  &&  (u <  k_right[it])}[0]
    }
>>>>>>> e36215d7f928cdf5d229fb91068c00a8c73570e8

    operator fun invoke(u: FloatArray): Array<FloatArray>{

<<<<<<< HEAD
    }

    operator fun invoke(u: Float): FloatArray {
        val index = bin(u, this.k)
        return deBoor(u, index, this.k, this.c, this.p)
=======
        val s = FloatArray(u.size)
        for (i in u.indices) {
            val index = findKnotIndex(u[i], this.k)
            val v  = (u[i] - this.k[index]) / (this.k[index] - this.k[index + 1])
            s[i] = evalBezier(v, this.b[indexgit])
>>>>>>> e36215d7f928cdf5d229fb91068c00a8c73570e8

    }
}
val l = linspace(0f,1f,4)

fun main(){

<<<<<<< HEAD
}

=======
    val spline = Spline(x, y, u)
    println(spline.findKnotIndex(0.5f, l))
}
>>>>>>> e36215d7f928cdf5d229fb91068c00a8c73570e8
