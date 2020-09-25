package com.e.minipaint

import java.util.*
import kotlin.math.pow

/*
Here we collect numpy-like methods for FloatArray.
 */
fun FloatArray.pow(p:Float) = FloatArray(this.size){ i -> this[i].pow(p)}
fun FloatArray.pow(n:Int) = FloatArray(this.size){ i -> this[i].pow(n)}
val Boolean.float get() = if (this) 1F else 0F //In order to do math with the truth-value of an expression.
// Might be better to use if-statement in the lambda
fun FloatArray.toString() = Arrays.toString(this)
fun IntArray.toString() =  Arrays.toString(this)
operator fun FloatArray.minus(fl: Float): FloatArray {return FloatArray(this.size, {i -> this[i] - fl})}
operator fun FloatArray.div(fl: Float): FloatArray {return FloatArray(this.size, {i -> this[i] / fl})}
operator fun FloatArray.times(a: FloatArray) = FloatArray(this.size){ i -> this[i] * a[i]}
operator fun FloatArray.times(a: Float) = FloatArray(this.size){ i -> this[i] * a }
operator fun Float.times(a: FloatArray) = FloatArray(a.size){ i -> a[i] * this }


fun norm(a: FloatArray, p: Float = 2.0F ): Float{
    return a.pow(p).sum().pow(1/p)
}

fun linspace(a: Float, b: Float, n: Int):FloatArray{
    return FloatArray(n){i -> a * (n - 1 - i ) / (n - 1) + b * i / (n - 1)}
}

fun cumsum(a: FloatArray): FloatArray{
    return FloatArray(a.size, {i -> a.slice(0..i+1).sum()})
}

fun bin(a: FloatArray, k : FloatArray): IntArray {
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
fun bin(a: Float, k : FloatArray): Int {
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

class A(a: FloatArray){
    val x: FloatArray = a
    override fun toString(): String {
        return Arrays.toString(this.x)

    }
    operator fun invoke(): Int{
        return x.size
    }
}

val a = A(floatArrayOf(1.0F, 20.0F))

fun printA(a: FloatArray){
    println(Arrays.toString(a))
}
fun printA(a: IntArray){
    println(Arrays.toString(a))
}
fun printA(a: List<FloatArray>){
    for (x in a) println("${Arrays.toString(x)},")
}

fun solveTridiagonal(A : List<FloatArray>,  b : FloatArray): FloatArray {
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
    var w : Float
    val n = b.size
    for (i in 1 .. n - 1) {
        w = a[i - 1]/b[i - 1]
        b[i] = b[i] - w * c[i - 1]
        d[i] = d[i] - w * d[i - 1]
    }
    val x = FloatArray(n)
    x[n-1] = d[n - 1]/b[n-1]
    for (i in n-2 downTo 0) {
        x[i] = (d[i] - c[i] * x[i + 1])/b[i]
    }
    return x
}