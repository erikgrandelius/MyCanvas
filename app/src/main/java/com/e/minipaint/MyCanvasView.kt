package com.e.minipaint

import android.content.Context
import android.graphics.Bitmap
import android.view.View
import androidx.core.content.res.ResourcesCompat
import android.graphics.Canvas
import android.graphics.Paint
import android.graphics.Path
import android.view.MotionEvent
import android.view.ViewConfiguration

import android.graphics.Matrix
import java.util.*


class MyCanvasView(context: Context): View(context) {
    private lateinit var extraCanvas: Canvas
    private lateinit var extraBitmap: Bitmap
    private val backgroundColor = ResourcesCompat.getColor(resources, R.color.colorBackground, null)
    private val drawColor = ResourcesCompat.getColor(resources, R.color.colorPaint, null)
    companion object {
        private const val STROKE_WIDTH = 6f // has to be float
    }
    private val paint = Paint().apply {
        color = drawColor
        // Smooths out edges of what is drawn without affecting shape.
        isAntiAlias = true
        // Dithering affects how colors with higher-precision than the device are down-sampled.
        isDither = true
        style = Paint.Style.STROKE // default: FILL
        strokeJoin = Paint.Join.ROUND // default: MITER
        strokeCap = Paint.Cap.ROUND // default: BUTT
        strokeWidth = STROKE_WIDTH // default: Hairline-width (really thin)
    }

    private var path = Path()
    private var motionTouchEventX = 0f
    private var motionTouchEventY = 0f
    private var currentX = 0f
    private var currentY = 0f

    var rawStrokeX = mutableListOf<Float>()
    var rawStrokeY = mutableListOf<Float>()



    private val touchTolerance = ViewConfiguration.get(context).scaledTouchSlop * 00.1

    override fun onSizeChanged(width: Int, height: Int, oldWidth: Int, oldHeight: Int) {
        super.onSizeChanged(width, height, oldWidth, oldHeight)
        if (::extraBitmap.isInitialized) {
            extraBitmap.recycle()
        }
        extraBitmap = Bitmap.createBitmap(width, height, Bitmap.Config.ARGB_8888)
        extraCanvas = Canvas(extraBitmap)
        extraCanvas.drawColor(backgroundColor)
    }

    override fun onDraw(canvas: Canvas) {
        super.onDraw(canvas)
        canvas.drawBitmap(extraBitmap, 0f, 0f, null)
    }

    override fun onTouchEvent(event: MotionEvent): Boolean {
        motionTouchEventX = event.x
        motionTouchEventY = event.y

        when (event.action) {
            MotionEvent.ACTION_DOWN -> touchStart()
            MotionEvent.ACTION_MOVE -> touchMove()
            MotionEvent.ACTION_UP -> touchUp()
        }
        return true
    }
    private fun touchStart() {
        path.reset()
        path.moveTo(motionTouchEventX, motionTouchEventY)
        println(motionTouchEventX)
        currentX = motionTouchEventX
        currentY = motionTouchEventY
    }

    private fun touchMove() {
        val dx = Math.abs(motionTouchEventX - currentX)
        val dy = Math.abs(motionTouchEventY - currentY)
        if (dx >= touchTolerance || dy >= touchTolerance) {
            path.lineTo(motionTouchEventX, motionTouchEventY)
            currentX = motionTouchEventX
            currentY = motionTouchEventY
            rawStrokeX.add(currentX)
            rawStrokeY.add(currentY)
            // Draw the path in the extra bitmap to cache it.
            extraCanvas.drawPath(path, paint)
        }
        invalidate()
    }

    private fun touchUp() {
        val model = Spline(rawStrokeX.toFloatArray(), rawStrokeY.toFloatArray(), rawStrokeX.toFloatArray())
        val a = model.b
        path.reset()
        val cubicPath = Path()
        cubicPath.moveTo(rawStrokeX[0], rawStrokeY[0])
        for (i in 0..rawStrokeX.size - 2){
            cubicPath.cubicTo(a[0][i], a[1][i], a[2][i], a[3][i], rawStrokeX[i +1], rawStrokeY[i + 1])
        }
        extraCanvas.drawPath(cubicPath, paint)
        rawStrokeX.clear()
        rawStrokeY.clear()
    }

}