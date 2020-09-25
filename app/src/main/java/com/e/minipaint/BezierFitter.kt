package com.e.minipaint


import java.lang.Math.sqrt
import java.util.*
import kotlin.math.pow

/*
Here we collect numpy-like methods for DoubleArray.
 */
//In order to do math with the truth-value of an expression. Might be better to use if-statement in the lambda
val Boolean.double get() = if (this) 1.0 else 0.0
fun DoubleArray.toString() = Arrays.toString(this)
fun IntArray.toString() =  Arrays.toString(this)
operator fun DoubleArray.times(a: DoubleArray) = DoubleArray(this.size){ i -> this[i] * a[i]}
operator fun DoubleArray.times(a: Double) = DoubleArray(this.size){ i -> this[i] * a }
operator fun Double.times(a: DoubleArray) = DoubleArray(a.size){ i -> a[i] * this }
fun DoubleArray.pow(p:Double) = DoubleArray(this.size){i -> this[i].pow(p)}
fun DoubleArray.pow(n:Int) = DoubleArray(this.size){i -> this[i].pow(n)}

fun norm(a: DoubleArray, p: Double = 2.0 ): Double{
    return a.pow(p).sum().pow(1/p)
}

fun linspace(a: Double, b: Double, n: Int):DoubleArray{
    return DoubleArray(n){i -> a * (n - 1 - i ) / (n - 1) + b * i / (n - 1)}
}

fun bin(a: DoubleArray, k : DoubleArray): IntArray {
    /*
    For values x and bin k, return Array of indices j such that
    k[J[j]]<= x[J[j]] <= k[J[j] + 1 ]
    @param x, array of values
    @param k, ordered array of bins

    returns: J, array of indices
     */
    val indices = IntArray(a.size)
    for (i in a.indices) {
        for (j in k.indices){
            if (a[i] < k[j]){
                indices[i] = j
                break
            }
            else{
                if(a[i] == k.last()){ indices[i] = k.size -1} else{indices[i] = k.size} }
        }

    }
    return indices
}
fun bin(a: Double, k : DoubleArray): Int {
    /*
    For values x and bin k, return Array of indices j such that
    k[J[j]]<= x[J[j]] < k[J[j] + 1 ]
    @param x, array of values
    @param k, ordered array of bins

    returns: J, array of indices
     */
    var index : Int = 0
    for (j in k.indices){
        if (a < k[j]){
            index = j
            break
            }
        else{
            if(a == k.last()){ index = k.size -1} else{index = k.size} }
    }

    return index
}



class A(a: DoubleArray){
    val x: DoubleArray = a
    override fun toString(): String {
        return Arrays.toString(this.x)

    }
    operator fun invoke(): Int{
        return x.size
    }
}

val a = A(doubleArrayOf(1.0, 20.0))

fun printA(a: DoubleArray){
    println(Arrays.toString(a))
}
fun printA(a: IntArray){
    println(Arrays.toString(a))
}
fun printA(a: List<DoubleArray>){
    for (x in a) println("${Arrays.toString(x)},")
}

val k = doubleArrayOf(0.0, 0.0, 0.0, 0.0, 1.0 ,2.0,3.0,4.0, 4.0, 4.0 ,4.0 )
val c_0 = doubleArrayOf(0.0, 0.0)
val c_1 = doubleArrayOf(1.0, 0.0)
val c_2 = doubleArrayOf(1.0, 1.0)
val c_3 = doubleArrayOf(1.0, 0.0)
val c_4 = doubleArrayOf(0.0, 0.0)
val c = mutableListOf<DoubleArray>(c_0, c_1, c_2, c_3,c_4)

val u = doubleArrayOf(0.0, 1.0, 2.0, 3.0, 4.0)


fun main(){
    val S = Spline(k,c, u, 3)
    println(S(2.0))
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
    fun solveTridiagonal(A : List<DoubleArray>,  b : DoubleArray): DoubleArray {
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
        val a = A[0]
        val b = A[1]
        val c = A[2]
        var w : Double
        val n = b.size
        for (i in 1 .. n - 1) {
            w = a[i - 1]/b[i - 1]
            b[i] = b[i] - w * c[i - 1]
            d[i] = d[i] - w * d[i - 1]
        }
        val x = DoubleArray(n)
        x[n-1] = d[n - 1]/b[n-1]
        for (i in n-2 downTo 0) {
            x[i] = (d[i] - c[i] * x[i + 1])/b[i]
        }
        return x
    }

    fun getBezierCoef(x: DoubleArray, y: DoubleArray): List<DoubleArray>{
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
        val C = listOf(DoubleArray(n-1, { i -> 1 + (i == n-1).double}),
            DoubleArray(n, {i -> 4 - 2 * (i == 0).double + 3 * (i == n).double}),
            DoubleArray(n- 1, {i -> 1.0})
        )
        val P_x = DoubleArray(n, { i -> 2 * (2 * x[i] + x[i + 1])})
        P_x[0] = x[0] + 2 * x[1]
        P_x[n-1] = 8 * x[n - 1] + x[n]
        val P_y = DoubleArray(n, {i -> 2 * (2 * y[i] + y[i + 1])})
        P_y[0] = y[0] + 2 * y[1]
        P_y[n-1] = 8 * y[n - 1] + y[n]
        val A_x = this.solveTridiagonal(C, P_x)
        val A_y = this.solveTridiagonal(C, P_y)
        val B_x = DoubleArray(n, { i -> if (i < n-1) 2 * x[i + 1] - A_x[i + 1] else (A_x[n - 1] + x[n]) / 2})
        val B_y = DoubleArray(n, { i -> if (i < n-1) 2 * y[i + 1] - A_y[i + 1] else (A_y[n - 1] + y[n]) / 2})

        return listOf(A_x, A_y, B_x, B_y)
    }
}

class Spline(val k: DoubleArray,val c: MutableList<DoubleArray>, val u : DoubleArray, val p: Int) {
    /* Class for implementing a cubic spline in B-spline form, see https://en.wikipedia.org/wiki/B-spline.
     A cubic spline is a piecewise polynomial, with  continuous second derivatives. The spline S is
     represented as a linear combination of cardinal B-spline basis functions,
     S(u) = sum_i=1^n c_i B_i(u).
     The parameter u is the normalized arc length, u in [0,1]. The knot vector k is not necessarily
     equidistantly spaced, i.e. k_i - k_i+1 may vary with i.
     */

    fun interpolate(u: DoubleArray, x: Array<DoubleArray>){
        return
    }
    fun deBoor(u : Double, index: Int,  k : DoubleArray, c: List<DoubleArray>, p : Int): DoubleArray{
        /*
        Evaluate S(x). Implementation of deBoors algorithm, see
        https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
        @param x : Double. Point where to evaluate S(x)
        @param index:  The knot interval where x belongs  k[i] <= x < k[i+1]
        @param k : List<Double>. Knot vector. Should be padded in end points
        @param c : List<DoubleArray>. Control point vector
         */

    println(index)
    val d = MutableList(p+2){j -> c[j + index -3 - p]}
    printA(d)
    var alpha: Double = 0.0
        for (r in 1 .. p) {
            for (j in p+1 downTo r ) {
                println("r: $r, j: $j, index: $index")
                alpha = (u - k[j + index - p]) / (k[j + 1 + index - r] - k[j + index - p])
                d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j]
            }
        }

        return d[p]

    }
    operator fun invoke(u: Double): DoubleArray{
        val index = bin(u,this.k)
        return deBoor(u, index, this.k, this.c, this.p)

    }


}