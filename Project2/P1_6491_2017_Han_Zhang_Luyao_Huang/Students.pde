 //<>// //<>// //<>// //<>//
final int point_number = 41;
pt [] pointSetCircleArcL = new pt[point_number];
pt [] pointSetCircleArcR = new pt[point_number];
pt [][] pointSetLowerArcs = new pt [4][point_number];
int ImageArraySize = 51;

boolean [][] texInArea = new boolean [ImageArraySize][ImageArraySize];
pt [] pointSet_MedialAxis = new pt[point_number * 4];


ArrayList<LocalCoorArea> CoordinateR = new ArrayList<LocalCoorArea>();
pt tangentPointL = new pt();
pt tangentPointR = new pt();

int record = 0;
pt curveStartPoint;

CIRCLE Circle_NL, Circle_NR, Circle_ML, Circle_MR;

//the nubmer of intersection of 4 arc is 3
LocalCoordinates [][] LowerArcsBotLC = new LocalCoordinates[4][point_number - 1];
LocalCoordinates lastPoint_LC;

ImageTexture imageR = new ImageTexture();
pt [][] imageLocalCoor = new pt [ImageArraySize][ImageArraySize];

// calculate the local coordinates once 
boolean recordFlag = true;

class LocalCoordinates
{
  float ori_x, ori_y;
  float x, y;
  LocalCoordinates() {
  };
  LocalCoordinates(float ori_x_, float ori_y_, float x_, float y_)
  {
    ori_x = ori_x_;
    ori_y = ori_y_;
    x = x_;
    y = y_;
  }
  void setLocalCoordinates(float ori_x_, float ori_y_, float x_, float y_)
  {
    ori_x = ori_x_;
    ori_y = ori_y_;
    x = x_;
    y = y_;
  }
}

LocalCoordinates Lerp_LC(LocalCoordinates a, LocalCoordinates b, float t)
{
  LocalCoordinates resultlc = new LocalCoordinates();

  resultlc.x = a.x + t * (b.x - a.x);
  resultlc.y = a.y + t * (b.y - a.y);
  resultlc.ori_x = a.ori_x + t * (b.ori_x - a.ori_x);
  resultlc.ori_y = a.ori_y + t * (b.ori_y - a.ori_y);

  return resultlc;
}

class LocalCoorArea
{
  LocalCoordinates startPoint;
  LocalCoordinates [] allPointsInSameX;
  LocalCoorArea()
  {
    allPointsInSameX = new LocalCoordinates[point_number];
  }
}




pt CalculatePointInMA(LocalCoordinates point, CIRCLE oriCircle, CIRCLE firstCircle, CIRCLE secondCircle)
{
  pt resultPoint = new pt();


  return resultPoint;
}


class SixArcArea

{
  CircleArc [] upperArcs;
  CircleArc [] lowerArcs;  
  int xaxis_resolution;
  int yaxis_resolution;
  SixArcArea() {
  };
  void Init()
  {
    pt S=P.G[0], E=P.G[1], L=P.G[2], R=P.G[3];
    pt M=P.G[4], N=P.G[5];

    curveStartPoint = S;

    xaxis_resolution = 40;
    yaxis_resolution = 40;

    upperArcs = new CircleArc [4];
    lowerArcs = new CircleArc [4];

    float s=d(S, L), e=d(E, R); 
    CIRCLE Cs = C(S, s), Ce = C(E, e); 

    upperArcs[0] = new CircleArc(Cs, pointSetCircleArcL[point_number / 2], L);
    upperArcs[1] = new CircleArc(Circle_NL, L, N);
    upperArcs[2] = new CircleArc(Circle_NR, N, tangentPointR);
    upperArcs[3] = new CircleArc(Ce, tangentPointR, pointSetCircleArcR[point_number / 2]);

    lowerArcs[0] = new CircleArc(Cs, pointSetCircleArcL[point_number / 2], tangentPointL);
    lowerArcs[1] = new CircleArc(Circle_ML, tangentPointL, M);
    lowerArcs[2] = new CircleArc(Circle_MR, M, R);
    lowerArcs[3] = new CircleArc(Ce, R, pointSetCircleArcR[point_number / 2]);



    //float deltaDistance = lengthSumOfLowerArcs / xaxis_resolution;

    //test
    for (int i = 0; i < 4; i++)
    {
      drawSmallArc(upperArcs[i].startPoint, upperArcs[i].endPoint, upperArcs[i].circle, 1);  
      drawSmallArc(pointSetLowerArcs[i], lowerArcs[i].startPoint, lowerArcs[i].endPoint, lowerArcs[i].circle, 1);
    }

    float lengthSumOfLowerArcs = 0;
    for (int i = 0; i < 4; i++)
    {
      lengthSumOfLowerArcs += lowerArcs[i].Length();
    }

    float startX = 0;

    for (int i = 0; i < 4; i++)
    {

      //test
      // float localdeltaxtest = (lowerArcs[i].Length() / lengthSumOfLowerArcs) / (point_number - 1);
      // println(localdeltaxtest * 40);

      for (int j = 0; j < point_number - 1; j++)
      {

        // LowerArcsBotLC[i * (point_number - 1) + j] = new LocalCoordinates();

        float localdeltax = (lowerArcs[i].Length() / lengthSumOfLowerArcs) / (point_number - 1);

        float local_x = startX + localdeltax * j;
        float local_y = -1;

        LowerArcsBotLC[i][j] =  new LocalCoordinates();
        LowerArcsBotLC[i][j].setLocalCoordinates(pointSetLowerArcs[i][j].x, pointSetLowerArcs[i][j].y, local_x, local_y);
        //test
        //  println(LowerArcsBotLC[i][j].ori_x,LowerArcsBotLC[i][j].ori_y);
      }

      startX += lowerArcs[i].Length() / lengthSumOfLowerArcs;
    }
    //last point;
    lastPoint_LC = new LocalCoordinates(pointSetCircleArcR[point_number / 2].x, pointSetCircleArcR[point_number / 2].y, 1, -1);

    //Special Condition, if the arc calculated is original circle arc(Cs,Ce), the point_number will change
    float localdeltax_s1 = (lowerArcs[0].Length() / lengthSumOfLowerArcs) / (point_number - 1);
    for (int i = point_number/ 2 + 1; i < point_number - 1; i++)
    {
      pt pointOnArcA = pointSetCircleArcL[i];

      //startPoint_x = 0

      float local_x = localdeltax_s1 * (i - point_number/ 2);
      float local_y = -1;

      LowerArcsBotLC[0][i - point_number/ 2 - 1] =  new LocalCoordinates();
      LowerArcsBotLC[0][i - point_number/ 2 - 1].setLocalCoordinates(pointOnArcA.x, pointOnArcA.y, local_x, local_y);
    }

    float localdeltax_s2 = (lowerArcs[3].Length() / lengthSumOfLowerArcs) / (point_number - 1);
    float startX_s2 = (lowerArcs[0].Length() + lowerArcs[1].Length() + lowerArcs[2].Length()) / lengthSumOfLowerArcs;
    for (int i = 1; i < (point_number - 1) / 2; i++)
    {
      pt pointOnArcA = pointSetCircleArcR[i];

      //startPoint_x = 0

      float local_x = startX_s2 + localdeltax_s2 * i;
      float local_y = -1;

      LowerArcsBotLC[3][i - 1] =  new LocalCoordinates();
      LowerArcsBotLC[3][i - 1].setLocalCoordinates(pointOnArcA.x, pointOnArcA.y, local_x, local_y);
    }





    //test
    /* beginShape();
     for(int i = 0;i < 4;i++)
     for(int j = 0;j < point_number -1;j++)
     {
     vertex(LowerArcsBotLC[i][j].ori_x,LowerArcsBotLC[i][j].ori_y);
     // vertex(pointSetLowerArcs[i][j].x,pointSetLowerArcs[i][j].y);
     // println(pointSetLowerArcs[i][j].x,pointSetLowerArcs[i][j].y);
     
     }
     endShape();*/
  }
}
// Only for inferior arc
class CircleArc
{
  CIRCLE circle;
  pt startPoint, endPoint;
  float startAngle, endAngle;

  CircleArc() {
  };
  CircleArc(CIRCLE circle_, pt spt, pt ept)
  {
    circle = circle_;
    startPoint = spt;
    endPoint = ept;

    startAngle = atan2((startPoint.y - circle.C.y)/circle.r, (startPoint.x - circle.C.x) /circle.r);
    endAngle = atan2((endPoint.y - circle.C.y) /circle.r, (endPoint.x - circle.C.x) /circle.r);
  }
  void setCricleArc(CIRCLE circle_, pt spt, pt ept)
  {
    circle = circle_;
    startPoint = spt;
    endPoint = ept;
    startAngle = atan2((startPoint.y - circle.C.y)/circle.r, (startPoint.x - circle.C.x) /circle.r);
    endAngle = atan2((endPoint.y - circle.C.y) /circle.r, (endPoint.x - circle.C.x) /circle.r);
  }

  float Length()
  {
    float angle = abs(endAngle - startAngle);
    angle = angle > PI? 2 * PI - angle: angle;
    float arcLength = angle * circle.r;
    return arcLength;
  }
  // only useful when point in on the circle of this arc
  boolean isPointOnArc(pt point)
  {
    vec vstart = new vec(startPoint.x - circle.C.x, startPoint.y - circle.C.y); 
    vec vend = new vec(endPoint.x - circle.C.x, endPoint.y - circle.C.y); 
    vec pvec = new vec(point.x - circle.C.x, point.y - circle.C.y); 

    float angleS = angle(pvec, vstart);
    float angleE = angle(pvec, vend);

    // println(angleS);
    // println(angleE);
    //  println(angleS * angleE <= 0 );
    if (angleS * angleE <= 0.00001 && abs(angleS) + abs(angleE) < PI)
      return true;
    else
      return false;
  }
}

//Same as from global to local, only the result is different
boolean inAreaR(pt globalpt)
{

  boolean foundflag = false;
  //look for x
  for (int i = 0; i < CoordinateR.size() - 1; i++)
  {
    for (int j = 0; j < point_number - 1; j++)
    {
      pt [] quad = new pt[4];
      quad[0] = new pt(CoordinateR.get(i).allPointsInSameX[j].ori_x, CoordinateR.get(i).allPointsInSameX[j].ori_y);
      quad[1] = new pt(CoordinateR.get(i + 1).allPointsInSameX[j].ori_x, CoordinateR.get(i  + 1).allPointsInSameX[j].ori_y);
      quad[2] = new pt(CoordinateR.get(i + 1).allPointsInSameX[j + 1].ori_x, CoordinateR.get(i + 1).allPointsInSameX[j + 1].ori_y);
      quad[3] = new pt(CoordinateR.get(i).allPointsInSameX[j + 1].ori_x, CoordinateR.get(i).allPointsInSameX[j + 1].ori_y);

      if (inQuadGlobal(globalpt, quad))
      {
        foundflag = true;
        break;
      }
    }
    if (foundflag)
      break;
  }

  if (foundflag)
  {
    return true;
  }
  return false;
}
//clockwise
boolean inQuadGlobal(pt testpoint, pt [] quadPoints)
{
  int intersectionCount = 0;

  for (int i = 0; i < 4; i++)
  {
    float t = (testpoint.x - quadPoints[i].x) / (quadPoints[(i + 1) % 4 ].x - quadPoints[i].x);
    float y = quadPoints[i].y + t * (quadPoints[(i + 1) % 4 ].y - quadPoints[i].y);
    if (t > -0.0001 && t < 1.0001 && y < testpoint.y)
      intersectionCount++;
  }

  if (intersectionCount % 2 == 0)
    return false;
  else
    return true;

  // float t2 = (testpoint.x - a.x) / (b.x - a.x);
  // float t3 = (testpoint.x - a.x) / (b.x - a.x);
  // float t4 = (testpoint.x - a.x) / (b.x - a.x);
}


pt LerpInQuadLocal(pt point, LocalCoordinates [] quad)
{
  pt resultpt = new pt();

  Line [] quadlines= new Line [4];
  for (int i = 0; i < 4; i++)
  {
    quadlines[i] = new Line();
    pt point_a = new pt(quad[i].ori_x, quad[i].ori_y);
    pt point_b = new pt(quad[(i + 1) % 4].ori_x, quad[(i + 1) % 4].ori_y);

    quadlines[i].figureoutLine(point_a, point_b);
    if (quadlines[i].onLine(point))
    {
      float t = (point.x - point_a.x) / (point_b.x - point_a.x);
      resultpt.x = quad[i].x + t * (quad[(i + 1) % 4].x - quad[i].x);
      resultpt.y = quad[i].y + t * (quad[(i + 1) % 4].y - quad[i].y);
      return resultpt;
    }
  }

  pt quadCorner = new pt(quad[0].ori_x, quad[0].ori_y);

  Line line = new Line();
  line.figureoutLine(quadCorner, point);

  pt intersectionPoint1 = intersection_LineAndLine(line, quadlines[1]);
  pt intersectionPoint2 = intersection_LineAndLine(line, quadlines[2]);

  // t1 and t2 are larger than 1
  float t1 = (point.x - quadCorner.x) / (intersectionPoint1.x - quadCorner.x);
  float t2 = (point.x- quadCorner.x) / (intersectionPoint2.x - quadCorner.x);

  pt intersectionPoint;
  float intersect_t;
  pt intersect_local = new pt();
  if (t1 > t2)
  {
    intersectionPoint = intersectionPoint1;
    intersect_t = t1;

    float t = (intersectionPoint.x - quad[1].ori_x) / (quad[2].ori_x - quad[1].ori_x);
    intersect_local.x =  quad[1].x + t * (quad[2].x - quad[1].x);
    intersect_local.y =  quad[1].y + t * (quad[2].y - quad[1].y);
  } else
  {
    intersectionPoint = intersectionPoint2;
    intersect_t = t2;

    float t = (intersectionPoint.x - quad[2].ori_x) / (quad[3].ori_x - quad[2].ori_x);
    intersect_local.x =  quad[2].x + t * (quad[3].x - quad[2].x);
    intersect_local.y =  quad[2].y + t * (quad[3].y - quad[2].y);
  }

  resultpt.x = quad[0].x + intersect_t * (intersect_local.x - quad[0].x);
  resultpt.y = quad[0].y + intersect_t * (intersect_local.y - quad[0].y);




  return resultpt;
}
pt FromGlobalToLocal(pt globalpt)
{
  pt localpt = new pt();
  boolean foundflag = false;
  //look for x
  for (int i = 0; i < CoordinateR.size() - 1; i++)
  {
    for (int j = 0; j < point_number - 1; j++)
    {
      pt [] quad = new pt[4];
      quad[0] = new pt(CoordinateR.get(i).allPointsInSameX[j].ori_x, CoordinateR.get(i).allPointsInSameX[j].ori_y);
      quad[1] = new pt(CoordinateR.get(i + 1).allPointsInSameX[j].ori_x, CoordinateR.get(i  + 1).allPointsInSameX[j].ori_y);
      quad[2] = new pt(CoordinateR.get(i + 1).allPointsInSameX[j + 1].ori_x, CoordinateR.get(i + 1).allPointsInSameX[j + 1].ori_y);
      quad[3] = new pt(CoordinateR.get(i).allPointsInSameX[j + 1].ori_x, CoordinateR.get(i).allPointsInSameX[j + 1].ori_y);

      if (inQuadGlobal(globalpt, quad))
      {
        LocalCoordinates [] quadLocal = new LocalCoordinates [4];
        quadLocal[0] = CoordinateR.get(i).allPointsInSameX[j];
        quadLocal[1] = CoordinateR.get(i + 1).allPointsInSameX[j];
        quadLocal[2] = CoordinateR.get(i + 1).allPointsInSameX[j + 1];
        quadLocal[3] = CoordinateR.get(i).allPointsInSameX[j + 1];

        localpt = LerpInQuadLocal(globalpt, quadLocal);
        foundflag = true;
        break;
      }
    }
    if (foundflag)
      break;
  }

  if (!foundflag)
  {
    localpt = new pt(-1, -1);
  }


  return localpt;
}


pt FromlocalToGlobal(pt localpt)
{
  pt globalpt = new pt();

  boolean foundflag = false;

  float delta_y = 2.0f / (point_number - 1);

  //look for x
  for (int i = 0; i < CoordinateR.size() - 1; i++)
  {
    float x1 = CoordinateR.get(i).startPoint.x;
    float x2 = CoordinateR.get(i + 1).startPoint.x;

    foundflag = true;
    //bug possible
    //found the point
    if (localpt.x >= x1 && localpt.x <= x2)
    {
      int yIndex = (int)((localpt.y + 1.0f) / delta_y);
      int yIndex2 = yIndex + 1;
      //upborder
      if (yIndex == point_number - 1)
      {
        yIndex2 = yIndex;
        yIndex--;
      }
      float tx = (localpt.x - x1) / (x2 - x1);
      float ty = (localpt.y - CoordinateR.get(i).allPointsInSameX[yIndex].y) / (CoordinateR.get(i).allPointsInSameX[yIndex2].y - CoordinateR.get(i).allPointsInSameX[yIndex].y);

      LocalCoordinates lc1 = Lerp_LC(CoordinateR.get(i).allPointsInSameX[yIndex], CoordinateR.get(i + 1).allPointsInSameX[yIndex], tx);
      LocalCoordinates lc2 = Lerp_LC(CoordinateR.get(i).allPointsInSameX[yIndex2], CoordinateR.get(i + 1).allPointsInSameX[yIndex2], tx);
      LocalCoordinates result_lc = Lerp_LC(lc1, lc2, ty);

      globalpt.x = result_lc.ori_x;
      globalpt.y = result_lc.ori_y;

      break;
    }
  }

  if (!foundflag)
  {
    globalpt = new pt(-1, -1);
  }

  return globalpt;
}

void CreateCoordinateR()
{
  pt S=P.G[0], E=P.G[1], L=P.G[2], R=P.G[3];
  pt M=P.G[4], N=P.G[5];

  SixArcArea areaR = new SixArcArea();
  areaR.Init();

  CoordinateR.clear();

  // divide the two cricles into 4 arcs
  //First arc
  TraversalExtend(pointSetCircleArcL, S, LowerArcsBotLC[0], 1);
  //Seoncd and third arc
  int result_i = CreateMedialAxis(areaR.lowerArcs[1], areaR.upperArcs[1], pointSetLowerArcs[1], 0, 0, LowerArcsBotLC[1]);
  //next bot arc

  if (result_i == point_number)
  {
    result_i = CreateMedialAxis(areaR.lowerArcs[2], areaR.upperArcs[1], pointSetLowerArcs[2], 0, 0, LowerArcsBotLC[2]);

    result_i = CreateMedialAxis(areaR.lowerArcs[2], areaR.upperArcs[2], pointSetLowerArcs[2], result_i, 0, LowerArcsBotLC[2]);
  } else
  {
    result_i = CreateMedialAxis(areaR.lowerArcs[1], areaR.upperArcs[2], pointSetLowerArcs[1], result_i, 0, LowerArcsBotLC[1]);

    result_i = CreateMedialAxis(areaR.lowerArcs[2], areaR.upperArcs[2], pointSetLowerArcs[2], 0, 0, LowerArcsBotLC[2]);
  }

  //lastarc
  TraversalExtend(pointSetCircleArcR, E, LowerArcsBotLC[3], 0);

  if (recordFlag)
  {
    recordFlag = false;
    imageR.Init();
    
    for (int i = 0; i < ImageArraySize; i++)
      for (int j = 0; j < ImageArraySize; j++)
      {
        if(inAreaR(imageR.coordinates[i][j].position))
        {
          texInArea[i][j] = true;
          imageLocalCoor[i][j] = new pt();
          imageLocalCoor[i][j] = FromGlobalToLocal(imageR.coordinates[i][j].position);
        }
        
      }
  }
  
   for (int i = 0; i < ImageArraySize; i++)
      for (int j = 0; j < ImageArraySize; j++)
      {
        if(texInArea[i][j])
        {
          imageR.coordinates[i][j].position = FromlocalToGlobal(imageLocalCoor[i][j]);
        }
      }


  imageR.drawImage(); //<>//

  //test
  //println(CoordinateR.size());
  /*beginShape();
   for(int j = 0; j < CoordinateR.size();j++)
   for(int i = 0; i < point_number ;i++)
   vertex(CoordinateR.get(j).allPointsInSameX[i].ori_x,CoordinateR.get(j).allPointsInSameX[i].ori_y);
   endShape();*/
  //test
  /* println("***************************************************");
   for(int j = 0; j < CoordinateR.size();j++)
   {
   for (int i = 0; i < point_number; i++)
   {
   float x = CoordinateR.get(j).allPointsInSameX[i].x;
   float y = CoordinateR.get(j).allPointsInSameX[i].y;
   print(x + "," + y + "     ");
   }
   println("");
   }*/


  /* line(pointSetCircleArcL[point_number / 2].x,pointSetCircleArcL[point_number / 2].y,L.x,L.y);
   line(pointSetCircleArcL[point_number / 2].x,pointSetCircleArcL[point_number / 2].y,tangentPointL.x,tangentPointL.y);
   line(pointSetCircleArcR[point_number / 2].x,pointSetCircleArcR[point_number / 2].y,R.x,R.y);
   line(pointSetCircleArcR[point_number / 2].x,pointSetCircleArcR[point_number / 2].y,tangentPointR.x,tangentPointR.y);*/
  /* pt testCircleCenter = new pt(100, 100);
   CIRCLE testc = new CIRCLE(testCircleCenter, 100);
   pt testspt = testc.PtOnCircle(PI / 3);
   pt testept = testc.PtOnCircle(-PI * 2 / 3);
   CircleArc testarc = new CircleArc(testc,testspt,testept);
   println(testarc.Length());*/
}

void drawOriCircleArc()
{
  pt S=P.G[0], E=P.G[1], L=P.G[2], R=P.G[3];
  pt M=P.G[4], N=P.G[5];

  CIRCLE Cs = C(S, d(S, L)), Ce = C(E, d(E, R));

  vec vec_SL, vec_STL, vec_ER, vec_ETR;
  // vec vec_CTL, vec_CL, vec_CTR, vec_CR;

  vec_SL = new vec(L.x - S.x, L.y - S.y);
  vec_STL = new vec(N.x - S.x, N.y - S.y);

  vec_ER = new vec(R.x - E.x, R.y - E.y);
  vec_ETR = new vec(M.x - E.x, M.y - E.y);


  if (angle(vec_SL, vec_STL) < 0 )
  {
    pointSetonCircle(pointSetCircleArcL, L, tangentPointL, Cs, false);
  } else
  {
    pointSetonCircle(pointSetCircleArcL, L, tangentPointL, Cs, true);
  }

  if (angle(vec_ER, vec_ETR) < 0 )
  {
    pointSetonCircle(pointSetCircleArcR, R, tangentPointR, Ce, false);
  } else
  {
    pointSetonCircle(pointSetCircleArcR, R, tangentPointR, Ce, true);
  }
  // pointSetCircleArcR[point_number - 1] = tangentPointR;
  // line(pointSetCircleArcR[point_number - 1].x,pointSetCircleArcR[point_number - 1].y,tangentPointR.x,tangentPointR.y);

  pt CircleCenterNL = crossPoint((L.y-S.y) / (L.x-S.x), L.y - L.x*(L.y-S.y) / (L.x-S.x), -(N.x-L.x)/(N.y-L.y), (N.x+L.x)*(N.x-L.x)/(2*(N.y-L.y))+(N.y+L.y)/2);
  float radiusNL = d(CircleCenterNL, L);
  Circle_NL = new CIRCLE(CircleCenterNL, radiusNL);

  pt CircleCenterMR = crossPoint((R.y-E.y) / (R.x-E.x), R.y - R.x*(R.y-E.y) / (R.x-E.x), -(M.x-R.x)/(M.y-R.y), (M.x+R.x)*(M.x-R.x)/(2*(M.y-R.y))+(M.y+R.y)/2);
  float radiusMR = d(CircleCenterMR, R);
  Circle_MR = new CIRCLE(CircleCenterMR, radiusMR);

  /*beginShape();
   for (int i = 0; i < point_number; i++)
   {
   vertex(pointSetCircleArcL[i].x, pointSetCircleArcL[i].y);
   }
   endShape();
   
   beginShape();
   for (int i = 0; i < point_number; i++)
   {
   vertex(pointSetCircleArcR[i].x, pointSetCircleArcR[i].y);
   }
   endShape();*/
}
void drawArc(pt A, pt S, pt B, float b)
{
  pt O = new pt(0, 0);
  float r = (b*b-d(B, S)*d(B, S))/(2*dot(V(B, S), U(R(V(A, S)))));
  O = P(S, (U(R(V(A, S))).scaleBy(r)));
  pt E = new pt(0, 0);
  E = P(O, U(R(V(O, B), -acos(r/d(O, B)))).scaleBy(r));
  drawCircleArcInHat(S, O, E);
}

void drawArcFromPoint(pt A, pt S, pt B, float b, int leftOrRight)
{
  pt O = new pt(0, 0);
  float r = (b*b-d(B, S)*d(B, S))/(2*dot(V(B, S), U(V(A, S))));
  O = P(S, (U(V(A, S)).scaleBy(r)));
  pt E = new pt(0, 0);
  E = P(O, U(R(V(O, B), -acos(r/d(O, B)))).scaleBy(r));
  // drawCircleArcInHat(S, O, E);

  pt CircleCenter;
  CircleCenter = crossPoint(-(S.x - O.x) / (S.y - O.y), S.y - S.x * (-(S.x - O.x) / (S.y - O.y)), -(E.x - O.x) / (E.y - O.y), E.y - E.x * (-(E.x - O.x) / (E.y - O.y)));

  if (leftOrRight == 0)
  {
    Circle_ML = new CIRCLE(CircleCenter, d(CircleCenter, S));
    tangentPointL = E;

    //  Circle_ML.drawCirc();
  } else if (leftOrRight == 1)
  {
    Circle_NR = new CIRCLE(CircleCenter, d(CircleCenter, S));
    tangentPointR = E;
    // Circle_NR.drawCirc();
  }
}


void pointSetonCircle(pt [] pointSet, pt point_a, pt point_b, CIRCLE circle, boolean clockwise)
{
  float startAngle = atan2((point_a.y - circle.C.y)/circle.r, (point_a.x - circle.C.x) /circle.r);
  float endAngle = atan2((point_b.y - circle.C.y) /circle.r, (point_b.x - circle.C.x) /circle.r);

  // println(pow((point_b.y - circle.C.y) /circle.r,2) + "   " + pow((point_b.x - circle.C.x) /circle.r,2));

  if (clockwise)
  {
    if (endAngle > startAngle)
      startAngle += 2 * PI;
  } else
  {
    if (endAngle < startAngle)
      endAngle += 2 * PI;
  }

  float theta = startAngle;
  float deltaAngle = (endAngle - startAngle)  / (point_number - 1);
  for (int i = 0; i < point_number; i++, theta += deltaAngle)
  {
    float x = circle.C.x + circle.r * cos(theta);
    float y = circle.C.y + circle.r * sin(theta);
    pointSet[i] = new pt();
    pointSet[i].x = x;
    pointSet[i].y = y;
  }
}

void drawSmallArc(pt [] pointSet, pt point_a, pt point_b, CIRCLE circle, int smallArc)
{
  pointSetonCircle(pointSet, point_a, point_b, circle, smallArc);

  beginShape();
  for (int i = 0; i < point_number; i++)
  {
    vertex(pointSet[i].x, pointSet[i].y);
  }
  endShape();
}

void drawSmallArc(pt point_a, pt point_b, CIRCLE circle, int smallArc)
{
  pt [] pointSet = new pt [point_number];

  pointSetonCircle(pointSet, point_a, point_b, circle, smallArc);

  beginShape();
  for (int i = 0; i < point_number; i++)
  {
    vertex(pointSet[i].x, pointSet[i].y);
  }
  endShape();
}

void pointSetonCircle(pt [] pointSet, pt point_a, pt point_b, CIRCLE circle, int smallArc)
{


  float startAngle = atan2((point_a.y - circle.C.y)/circle.r, (point_a.x - circle.C.x) /circle.r);
  float endAngle = atan2((point_b.y - circle.C.y) /circle.r, (point_b.x - circle.C.x) /circle.r);

  // println(pow((point_b.y - circle.C.y) /circle.r,2) + "   " + pow((point_b.x - circle.C.x) /circle.r,2));

  /*if (endAngle < startAngle)
   {
   float temp = endAngle;
   endAngle = startAngle;
   startAngle = temp;
   }*/

  if (smallArc == 1)
  {
    if (abs(endAngle - startAngle) > PI )
    {
      endAngle = endAngle < 0? 2 * PI + endAngle:endAngle;
      startAngle = startAngle < 0? 2 * PI + startAngle:startAngle;
      /* float temp = endAngle;
       endAngle = startAngle + 2 * PI;
       startAngle = temp;*/
    }
  } else
  {
    // Won't be used in this project

    /*if (abs(endAngle - startAngle) < PI )
     {
     endAngle = endAngle > 0?endAngle - 2 * PI:endAngle;
                                  /*float temp = startAngle;
     startAngle = endAngle;
     endAngle = temp + 2 * PI;
     }*/
  }

  float theta = startAngle;
  float deltaAngle = (endAngle - startAngle)  / (point_number - 1);
  for (int i = 0; i < point_number; i++, theta += deltaAngle)
  {
    float x = circle.C.x + circle.r * cos(theta);
    float y = circle.C.y + circle.r * sin(theta);
    pointSet[i] = new pt();
    pointSet[i].x = x;
    pointSet[i].y = y;
  }
}

pt crossPoint(float k1, float m1, float k2, float m2)
{
  float x = (m2-m1)/(k1-k2);
  float y = k1*x+m1;
  pt cp = new pt(x, y);
  return cp;
}