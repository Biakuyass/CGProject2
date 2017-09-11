

class Matrix
{
  float [][] value;
  Matrix()
  {
    value = new float[3][3];
    value[0][0] = 1;
    value[0][1] = 0;
    value[0][2] = 0;
    value[1][0] = 0;
    value[1][1] = 1;
    value[1][2] = 0;
    value[2][0] = 0;
    value[2][1] = 0;
    value[2][2] = 1;
  }
  void setMatrix(Matrix b)
  {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
        value[i][j] = b.value[i][j];
      }
  }
  void setRotation(float theta)
  {
    value[0][0] = cos(theta);
    value[0][1] = sin(theta);
    value[0][2] = 0;
    value[1][0] = -sin(theta);
    value[1][1] = cos(theta);
    value[1][2] = 0;
    value[2][0] = 0;
    value[2][1] = 0;
    value[2][2] = 1;
  }

  void setRotation(float costheta, float sintheta)
  {
    value[0][0] = costheta;
    value[0][1] = sintheta;
    value[0][2] = 0;
    value[1][0] = -sintheta;
    value[1][1] = costheta;
    value[1][2] = 0;
    value[2][0] = 0;
    value[2][1] = 0;
    value[2][2] = 1;
  }


  void setTranslation(float x, float y)
  {
    value[0][0] = 1;
    value[0][1] = 0;
    value[0][2] = 0;
    value[1][0] = 0;
    value[1][1] = 1;
    value[1][2] = 0;
    value[2][0] = x;
    value[2][1] = y;
    value[2][2] = 1;
  }
}

Matrix Apply(Matrix mat1, Matrix mat2)
{
  Matrix resultMatrix = new Matrix();
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      resultMatrix.value[j][i] = mat1.value[j][0] * mat2.value[0][i] + mat1.value[j][1] * mat2.value[1][i] + mat1.value[j][2] * mat2.value[2][i];
    }
  return resultMatrix;
}

pt Apply(pt point, Matrix mat)
{
  pt resultpt = new pt();
  resultpt.x = point.x * mat.value[0][0] + point.y * mat.value[1][0] + mat.value[2][0];
  resultpt.y = point.x * mat.value[0][1] + point.y * mat.value[1][1] + mat.value[2][1];
  float scale = point.x * mat.value[0][2] + point.y * mat.value[1][2] + mat.value[2][2];

  resultpt.x /= scale;
  resultpt.y /= scale;

  return resultpt;
}

vec Apply(vec vector, Matrix mat)
{
  vec resultvec = new vec();

  resultvec.x = vector.x * mat.value[0][0] + vector.y * mat.value[1][0];
  resultvec.y = vector.x * mat.value[0][1] + vector.y * mat.value[1][1];
  float scale = vector.x * mat.value[0][2] + vector.y * mat.value[1][2] + mat.value[2][2];

  resultvec.x /= scale;
  resultvec.y /= scale;

  return resultvec;
}

void PrintMatrix(Matrix mat)
{
  for (int i = 0; i < 3; i++)
  {
    println(mat.value[i][0] + "  " + mat.value[i][1] + "  " + mat.value[i][2]);
  }
}

class Hyperbola
{
  float square_a, square_b, square_c;
  boolean isright; // decide left or right curve of Hyperbola;
  Hyperbola() {
  };
  Hyperbola(float square_a_, float square_b_, float square_c_) {
    square_a = square_a_; 
    square_b = square_b_; 
    square_c = square_c_;
  }

  void setHyperbola(float square_a_, float square_b_, float square_c_) {
    square_a = square_a_; 
    square_b = square_b_; 
    square_c = square_c_;
  }

  void pointSetInHyperbola(pt [] pointSet, pt point_a, pt point_b)
  {

    float a = sqrt(square_a);
    float b = sqrt(square_b);

    float startAngle = atan2((point_a.y * a) / (b * point_a.x), a / point_a.x);
    float endAngle = atan2((point_b.y * a) / (b * point_b.x), a / point_b.x);


    if (endAngle < startAngle)
    {
      float temp = endAngle;
      endAngle = startAngle;
      startAngle = temp;
    }

    if (endAngle - startAngle > PI)
    {
      float temp = endAngle;
      endAngle = startAngle + 2 * PI;
      startAngle = temp;
    }

    if (startAngle < PI / 2.0f && startAngle > -PI / 2.0f)
    {
      isright = true;
    } else
    {
      isright = false;
    }

    float theta = startAngle;
    float deltaAngle = (endAngle - startAngle)  / (point_number - 1);
    for (int i = 0; i < point_number; i++, theta += deltaAngle)
    {
      pointSet[i] = new pt();
      float x = sqrt(square_a) / cos(theta);
      float y = sqrt(square_b) * tan(theta);
      pointSet[i].x = x;
      pointSet[i].y = y;
    }
  }
}

class Ellipse
{
  float square_a, square_b, square_c;
  Ellipse() {
  };
  Ellipse(float square_a_, float square_b_, float square_c_) {
    square_a = square_a_; 
    square_b = square_b_; 
    square_c = square_c_;
  }

  boolean pointInArc(pt point_a, pt point_b, pt testPoint)
  {
    float a = sqrt(square_a);
    float b = sqrt(square_b);

    //vec vector_CLS = new vec(changed_S.x ,changed_S.y );
    // vec vector_CLE = new vec(changed_E.x ,changed_E.y );
    float startAngle = atan2(point_a.y / b, point_a.x / a);
    float endAngle = atan2(point_b.y / b, point_b.x / a);

    if (endAngle < startAngle)
    {
      float temp = endAngle;
      endAngle = startAngle;
      startAngle = temp;
    }

    if (endAngle - startAngle > PI)
    {
      float temp = endAngle;
      endAngle = startAngle + 2 * PI;
      startAngle = temp;
    }

    float testAngle = atan2(testPoint.y / b, testPoint.x / a);
    if (testAngle < endAngle && testAngle > startAngle )
      return true;
    else if (testAngle + 2 * PI < endAngle && testAngle + 2 * PI > startAngle )
      return true;
    else if (testAngle - 2 * PI < endAngle && testAngle + 2 * PI > startAngle )
      return true;

    return false;
  }

  void pointSetInEllipse(pt [] pointSet, pt point_a, pt point_b)
  {
    float a = sqrt(square_a);
    float b = sqrt(square_b);

    //vec vector_CLS = new vec(changed_S.x ,changed_S.y );
    // vec vector_CLE = new vec(changed_E.x ,changed_E.y );
    float startAngle = atan2(point_a.y / b, point_a.x / a);
    float endAngle = atan2(point_b.y / b, point_b.x / a);

    if (endAngle < startAngle)
    {
      float temp = endAngle;
      endAngle = startAngle;
      startAngle = temp;
    }

    if (endAngle - startAngle > PI)
    {
      float temp = endAngle;
      endAngle = startAngle + 2 * PI;
      startAngle = temp;
    }



    float theta = startAngle;
    float deltaAngle = (endAngle - startAngle)  / (point_number - 1);
    for (int i = 0; i < point_number; i++, theta += deltaAngle)
    {
      pointSet[i] = new pt();
      float x = sqrt(square_a) * cos(theta);
      float y = sqrt(square_b) * sin(theta);
      pointSet[i].x = x;
      pointSet[i].y = y;
    }
  }
}

class Line
{
  float k, b;
  boolean isperpendicular;
  float x_value; // x = x_value when isperpendicular is true
  Line()
  {
    k = 1;
    b = 1;
  };
  Line(float k_, float b_)
  {
    k = k_;
    b = b_;
  }

  void figureoutLine(pt point_a, pt point_b)
  {
    if (abs(point_b.x - point_a.x) < 0.0001)
    {
      isperpendicular = true;
      x_value = point_b.x;
    } else
    {
      isperpendicular = false;
      k = (point_b.y - point_a.y) / (point_b.x - point_a.x);
      b = point_a.y - k * point_a.x;
    }
  }
  
  boolean onLine(pt testPoint)
  {
    float resulty = testPoint.x * k + b;
    if(abs(resulty - testPoint.y) < 0.001f)
    return true;
    else
    return false;
    
  }
}
pt intersection_LineAndLine(Line line_a, Line line_b)
{
  pt resultpt = new pt();
  resultpt.x = (line_b.b - line_a.b) / (line_a.k - line_b.k);
  resultpt.y = line_a.k * resultpt.x + line_a.b;

  return resultpt;
}
//this method may cause a bug in special condition
//only for minor arc
pt intersection_LineAndArc(Line line, pt CircleCenter, float radius, pt line_startPoint)
{
  pt resultpt = new pt();

  float a = 1 + line.k * line.k;
  float b = 2 * line.k * line.b - 2 * CircleCenter.x - 2 * CircleCenter.y * line.k ;
  float c = line.b * line.b + CircleCenter.x * CircleCenter.x + CircleCenter.y * CircleCenter.y - 2 * CircleCenter.y * line.k - radius * radius;

  if (b * b - 4 * a * c < 0)
    println("Error! LineAndArc");
  float x1 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
  float y1 = line.k * x1 + line.b;
  pt pt1 = new pt(x1, y1);


  float x2 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
  float y2 = line.k * x2 + line.b;

  pt pt2 = new pt(x2, y2);

  resultpt = d(pt1, line_startPoint) < d(pt2, line_startPoint)? pt1 : pt2;

  return resultpt;
}

//this method may cause a bug in special condition
//only for minor ellipse arc
pt intersection_LineAndEllipse(Line line, Ellipse ellipse, pt line_startPoint)
{
  pt resultpt = new pt();

  float ellipse_a = sqrt(ellipse.square_a);
  float ellipse_b = sqrt(ellipse.square_b);

  float a = ellipse_b * ellipse_b + ellipse_a * ellipse_a * line.k * line.k;
  float b = 2 * ellipse_a * ellipse_a * line.k * line.b ;
  float c = ellipse_a * ellipse_a * line.b * line.b - ellipse_a * ellipse_a * ellipse_b * ellipse_b;

  if (b * b - 4 * a * c < 0)
    println("Error! bb Ellipse");

  float x1 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
  float y1 = line.k * x1 + line.b;
  pt pt1 = new pt(x1, y1);


  float x2 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
  float y2 = line.k * x2 + line.b;

  pt pt2 = new pt(x2, y2);

  resultpt = d(pt1, line_startPoint) < d(pt2, line_startPoint)? pt1 : pt2;

  return resultpt;
}

pt intersection_LineAndHyperbola(Line line_, Hyperbola hyperbola_, pt line_startPoint)
{


  pt resultpt = new pt();

  /*if(line_.isperpendicular)
   {
   resultpt.x = line_.x_value;
   resultpt.y = direction * sqrt(hyperbola_.square_b * ((line_.x_value * line_.x_value) / hyperbola_.square_a - 1));
   
   return resultpt;
   }*/


  float temp_a = hyperbola_.square_b - hyperbola_.square_a * line_.k * line_.k;
  float temp_b = -2 * hyperbola_.square_a * line_.k * line_.b;
  float temp_c = -1 * hyperbola_.square_b * hyperbola_.square_a - hyperbola_.square_a * line_.b * line_.b;


  float x1 = (-1 * temp_b + sqrt(sq(temp_b) - 4 * temp_a * temp_c) ) / (2 * temp_a);
  float y1 = line_.k * x1 + line_.b;
  pt point1 = new pt(x1, y1);

  float x2 = (-1 * temp_b - sqrt(sq(temp_b) - 4 * temp_a * temp_c) ) / (2 * temp_a);
  float y2 =line_.k * x2 + line_.b;
  pt point2 = new pt(x2, y2);

  if (hyperbola_.isright == true)
  {
    if (point1.x < 0.01)
      resultpt = point2;

    if (point2.x < 0.01)
      resultpt = point1;

    if (point1.x < 0.01 && point2.x < 0.01)
    {
      println("Error LineAndHyperbola");
      return new pt();
    }

    if (point1.x > -0.01 && point2.x > -0.01)
    {
      resultpt = d(point1, line_startPoint) < d(point2, line_startPoint) ? point1 : point2;
    }
  } else
  {
    if (point1.x > -0.01)
      resultpt = point2;

    if (point2.x > -0.01)
      resultpt = point1;

    if (point1.x > -0.01 && point2.x > -0.01)
    {
      println("Error LineAndHyperbola");
      return new pt();
    }

    if (point1.x < -0.01 && point2.x < -0.01)
    {
      resultpt = d(point1, line_startPoint) < d(point2, line_startPoint) ? point1 : point2;
    }
  }
  // resultpt.x = 


  // resultpt.y = line_.k * resultpt.x + line_.b;

  return resultpt;
}

Line figureoutTangentLine(pt circleCenter, pt point_)
{
  Line resultLine;
  resultLine = new Line();


  Line tempLine = new Line();
  tempLine.figureoutLine(circleCenter, point_);

  if (tempLine.isperpendicular)
  {
    resultLine.k = 0;
    resultLine.b = point_.y;

    return resultLine;
  }

  if (abs(tempLine.k) < 0.0001)
  {

    resultLine.isperpendicular = true;
    resultLine.x_value = point_.x;

    return resultLine;
  }


  float k = -1 / tempLine.k;
  float b = point_.y - k * point_.x;

  resultLine = new Line(k, b);

  return resultLine;
}

// intersection of two tangent line of a circle
pt tangentLineIntersection(pt circleCenter, pt point_a, pt point_b) 
{
  pt resultpt;

  Line tangentLineL, tangentLineR;

  tangentLineL = figureoutTangentLine(circleCenter, point_a);
  tangentLineR = figureoutTangentLine(circleCenter, point_b);

  resultpt = intersection_LineAndLine(tangentLineL, tangentLineR);

  return resultpt;
}


//clockwise
vec Rotate(vec vector, float theta)
{
  vec resultvec = new vec();
  Matrix matrix_T = new Matrix();
  matrix_T.setRotation(theta);
  resultvec = Apply(vector, matrix_T);
  return resultvec;
}

Matrix figureoutMatrix(float translate_x, float translate_y, float theta, float translate_x2, float translate_y2)
{
  Matrix translationMat = new Matrix();
  translationMat.setTranslation(translate_x, translate_y);

  Matrix rotationMat = new Matrix();
  rotationMat.setRotation(theta);

  Matrix translationMatAxis = new Matrix();
  translationMatAxis.setTranslation(translate_x2, translate_y2);

  Matrix tempMat = new Matrix();
  tempMat.setMatrix(Apply(translationMat, rotationMat));

  Matrix finalMat = new Matrix();
  finalMat.setMatrix(Apply(tempMat, translationMatAxis));

  return finalMat;
}

void LinePoint(pt point_a, pt point_b)
{
  line(point_a.x, point_a.y, point_b.x, point_b.y);
}

void TraversalExtend(pt [] pointSet, pt CircleCenter, LocalCoordinates [] localCoor, int leftorRight)
{
  pt line_endPoint = pointSet[point_number / 2];
  int division = (point_number - 1) / 2 - 1;
  float delta_t = 1.0f / division;

  for (int i = 1; i < (point_number - 1) / 2; i++)
  {
    pt pointOnMA = L(CircleCenter, line_endPoint, delta_t * i);
    pt pointOnArcA = pointSet[point_number - 1 - i];
    pt pointOnArcB = pointSet[i];

     
    //Calculate all coordiates when x = ....;
    LocalCoorArea oneCol = new LocalCoorArea();
    if (leftorRight == 1)
    {
      oneCol.startPoint = localCoor[(point_number - 1) / 2 - 1 - i];
      pt [] tempPointSet = new pt[point_number];
      GetGlobalYInHat(tempPointSet, pointOnArcA, pointOnMA, pointOnArcB);

      CalculateLocalY(oneCol.allPointsInSameX, tempPointSet, oneCol.startPoint);

      //line(oneCol.startPoint.ori_x,oneCol.startPoint.ori_y,pointOnMA.x,pointOnMA.y);

      CoordinateR.add(0,oneCol);
      
    } else
    {

      oneCol.startPoint = localCoor[i - 1];
      
      pt [] tempPointSet = new pt[point_number];
      GetGlobalYInHat(tempPointSet, pointOnArcB, pointOnMA, pointOnArcA);
      CalculateLocalY(oneCol.allPointsInSameX, tempPointSet, oneCol.startPoint);

      CoordinateR.add(oneCol);
    }




    //  LinePoint(pointOnArcA, pointOnMA);
    //   LinePoint(pointOnArcB, pointOnMA);
  }
}
//changed version of drawCircleArcInHat
void GetGlobalYInHat(pt[] pointSet, pt PA, pt B, pt PC) // draws circular arc from PA to PB that starts tangent to B-PA and ends tangent to PC-B
{
  float e = (d(B, PC)+d(B, PA))/2;
  pt A = P(B, e, U(B, PA));
  pt C = P(B, e, U(B, PC));
  vec BA = V(B, A), BC = V(B, C);
  float d = dot(BC, BC) / dot(BC, W(BA, BC));
  pt X = P(B, d, W(BA, BC));
  float r=abs(det(V(B, X), U(BA)));
  vec XA = V(X, A), XC = V(X, C); 

  float a = angle(XA, XC), da=a/(point_number - 1);
  //int i = 0;
  float w = 0;

  for (int i = 0; i < point_number; i++)
  {
    w = i * da;
    pointSet[i] = P(X, R(XA, w));
  }

  //  println(i);

  //beginShape(); 

  // endShape();
}  

void CalculateLocalY(LocalCoordinates [] localCoor, pt[] pointSet, LocalCoordinates startPoint)
{

  int division = point_number - 1;
  float delta_t = 2.0f / division;

  for (int i = 0; i < point_number; i++)
  {
    localCoor[i] = new LocalCoordinates();
    localCoor[i].ori_x = pointSet[i].x;
    localCoor[i].ori_y = pointSet[i].y;
    localCoor[i].x = startPoint.x;
    localCoor[i].y = startPoint.y + i * delta_t;
  }
}


int CreateMedialAxis(CircleArc arc_a, CircleArc arc_b, pt [] bottomPoints, int start_i, int test, LocalCoordinates [] localCoor)
{
  pt circleCenterL = arc_a.circle.C;
  pt circleCenterR = arc_b.circle.C;
  int divisionNumber = point_number;
  float radiusL = arc_a.circle.r;
  float radiusR = arc_b.circle.r;

  boolean contineue_flag = true;

  // pt S=P.G[0], E=P.G[1], L=P.G[2], R=P.G[3];
  // pt M=P.G[4], N=P.G[5];

  vec vector_LR = new vec(circleCenterR.x - circleCenterL.x, circleCenterR.y - circleCenterL.y);
  float theta = angle(vector_LR);

  Matrix finalMat = figureoutMatrix(-circleCenterL.x, -circleCenterL.y, -theta, -d(circleCenterL, circleCenterR)/2.0f, 0);

  Matrix finalMat_inv = figureoutMatrix(d(circleCenterL, circleCenterR)/2.0f, 0, theta, circleCenterL.x, circleCenterL.y);

  pt changed_circleCenterL, changed_circleCenterR, changed_S;


  changed_circleCenterL = Apply(circleCenterL, finalMat);
  changed_circleCenterR = Apply(circleCenterR, finalMat);
  changed_S = Apply(curveStartPoint, finalMat);


  float square_a, square_b, square_c;

  square_c = d2(changed_circleCenterL, changed_circleCenterR) * 0.25f;

  pt [] pointset = new pt [divisionNumber];
  //pt LastMAPoint = curveStartPoint;
  if ( abs(( d(curveStartPoint, circleCenterL) + d(curveStartPoint, circleCenterR) ) - (radiusL + radiusR)) < 0.01 )
  {

    square_a = pow((d(curveStartPoint, circleCenterL) + d(curveStartPoint, circleCenterR)) / 2.0f, 2);
    square_b = square_a - square_c;

    Ellipse temp_ellipse = new Ellipse(square_a, square_b, square_c); 

    int i = start_i;

    while (contineue_flag == true)
    {

      pt point = new pt(bottomPoints[i].x, bottomPoints[i].y);
      pt changed_Point = Apply(point, finalMat);

      Line line = new Line();
      line.figureoutLine(changed_Point, changed_circleCenterL);

      pt changed_pointOnMA = intersection_LineAndEllipse(line, temp_ellipse, changed_Point);

      Line line_b = new Line();
      line_b.figureoutLine(changed_pointOnMA, changed_circleCenterR);

      pt changed_pointOnArcB = intersection_LineAndArc(line_b, changed_circleCenterR, radiusR, changed_pointOnMA);




      pt pointOnArcB = Apply(changed_pointOnArcB, finalMat_inv);
      pt pointOnMA = Apply(changed_pointOnMA, finalMat_inv);



      if (!arc_b.isPointOnArc(pointOnArcB))
      {
        // line(pointOnMA.x,pointOnMA.y,point.x,point.y);
        //line(pointOnMA.x,pointOnMA.y,pointOnArcB.x,pointOnArcB.y);
        //println(i);
        contineue_flag = false;

        if (i > 0)
        {
          pt point2 = arc_b.endPoint;
          pt changed_Point2 = Apply(point2, finalMat);

          Line line2 = new Line();
          line2.figureoutLine(changed_Point2, changed_circleCenterR);

          pt changed_pointOnMA2 = intersection_LineAndEllipse(line2, temp_ellipse, changed_Point2);

          pt pointOnMA2 = Apply(changed_pointOnMA2, finalMat_inv);

          curveStartPoint = pointOnMA2;
        }
      }

      // println(i);
      if (contineue_flag)
      {
        if (i < point_number -1)
        {
          pointset[i] = pointOnMA;

          //Calculate all coordiates when x = ....;
          LocalCoorArea oneCol = new LocalCoorArea();
          oneCol.startPoint = localCoor[i];

          pt [] tempPointSet = new pt[point_number];
          GetGlobalYInHat(tempPointSet, point, pointOnMA, pointOnArcB);
          CalculateLocalY(oneCol.allPointsInSameX, tempPointSet, oneCol.startPoint);

          CoordinateR.add(oneCol);
        } else
          curveStartPoint = pointOnMA;

        i++;
        if (test == 1 && i < point_number -1 )
        {
          line(point.x, point.y, pointOnMA.x, pointOnMA.y);
          line(pointOnArcB.x, pointOnArcB.y, pointOnMA.x, pointOnMA.y);
        }


        if (i == point_number)
        {
          contineue_flag = false;
        }
      }
    }
    // println(i);
    beginShape();
    for (int j = start_i; j < i - 1; j++)
      vertex(pointset[j].x, pointset[j].y);
    endShape();

    return i;
    /*temp_ellipse.pointSetInEllipse(pointset, changed_S, changed_E);
     
     println("ellipse");
     
     for (int i = 0; i < divisionNumber; i++)
     {
     pointset[i] = Apply(pointset[i], finalMat_inv);
     }
     beginShape();
     
     for (int i = 0; i < divisionNumber; i++)
     {
     vertex(pointset[i].x, pointset[i].y);
     }
     endShape();*/
  } else if ( abs(abs(d(curveStartPoint, circleCenterL) - d(curveStartPoint, circleCenterR)) - abs(radiusL - radiusR)) < 0.01)
  {

    square_a = pow((d(curveStartPoint, circleCenterL) - d(curveStartPoint, circleCenterR)) / 2.0f, 2);
    square_b = square_c - square_a;

    float a = sqrt(square_a);
    float b = sqrt(square_b);


    Hyperbola temp_hyperbola = new Hyperbola(square_a, square_b, square_c); 

    float startAngle = atan2((changed_S.y * a) / (b * changed_S.x), a / changed_S.x);

    if (startAngle < PI / 2.0f && startAngle > -PI / 2.0f)
    {
      temp_hyperbola.isright = true;
    } else
    {
      temp_hyperbola.isright = false;
    }


    int i = start_i;

    while (contineue_flag == true)
    {

      pt point = new pt(bottomPoints[i].x, bottomPoints[i].y);
      pt changed_Point = Apply(point, finalMat);

      Line line = new Line();
      line.figureoutLine(changed_Point, changed_circleCenterL);

      pt changed_pointOnMA = intersection_LineAndHyperbola(line, temp_hyperbola, changed_Point);

      Line line_b = new Line();
      line_b.figureoutLine(changed_pointOnMA, changed_circleCenterR);

      pt changed_pointOnArcB = intersection_LineAndArc(line_b, changed_circleCenterR, radiusR, changed_pointOnMA);

      pt pointOnArcB = Apply(changed_pointOnArcB, finalMat_inv);
      pt pointOnMA = Apply(changed_pointOnMA, finalMat_inv);



      if (!arc_b.isPointOnArc(pointOnArcB))
      {
        // line(pointOnMA.x,pointOnMA.y,point.x,point.y);
        //line(pointOnMA.x,pointOnMA.y,pointOnArcB.x,pointOnArcB.y);
        //println(i);

        contineue_flag = false;

        if (i > 0)
        {
          pt point2 = arc_b.endPoint;
          pt changed_Point2 = Apply(point2, finalMat);

          Line line2 = new Line();
          line2.figureoutLine(changed_Point2, changed_circleCenterR);

          pt changed_pointOnMA2 = intersection_LineAndHyperbola(line2, temp_hyperbola, changed_Point2);

          pt pointOnMA2 = Apply(changed_pointOnMA2, finalMat_inv);

          curveStartPoint = pointOnMA2;
        }
      }

      // println(i);
      if (contineue_flag)
      {
        if (i < point_number -1)
        {

          pointset[i] = pointOnMA;

          //Calculate all coordiates when x = ....;
          LocalCoorArea oneCol = new LocalCoorArea();
          oneCol.startPoint = localCoor[i];

          pt [] tempPointSet = new pt[point_number];
          GetGlobalYInHat(tempPointSet, point, pointOnMA, pointOnArcB);
          CalculateLocalY(oneCol.allPointsInSameX, tempPointSet, oneCol.startPoint);

          CoordinateR.add(oneCol);
        } else
          curveStartPoint = pointOnMA;

        i++;
        if (test == 1 && i < point_number -1)
        {
          line(point.x, point.y, pointOnMA.x, pointOnMA.y);
          line(pointOnArcB.x, pointOnArcB.y, pointOnMA.x, pointOnMA.y);
        }

        //LastMAPoint = pointOnMA;

        if (i == point_number)
        {
          contineue_flag = false;
        }
      }
    }
    //println(i);
    beginShape();
    for (int j = start_i; j < i - 1; j++)
      vertex(pointset[j].x, pointset[j].y);
    endShape();

    return i;


    /*temp_hyperbola.pointSetInHyperbola(pointset, changed_S, changed_E);
     
     for (int i = 0; i < divisionNumber; i++)
     {
     pointset[i] = Apply(pointset[i], finalMat_inv);
     }
     
     
     beginShape();
     
     for (int i = 0; i < divisionNumber; i++)
     {
     vertex(pointset[i].x, pointset[i].y);
     }
     endShape();
     
     //circle part
     line(pointSetCircleArcL[point_number / 2].x, pointSetCircleArcL[point_number / 2].y, S.x, S.y);
     line(pointSetCircleArcR[point_number / 2].x, pointSetCircleArcR[point_number / 2].y, E.x, E.y);
     println("hyperbola0");*/
  } else
  {
    println("Special Condition");
    return start_i;
  }
}