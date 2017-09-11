//int ImageArraySize = 51;
class TextureVertex
{
  pt position;
  pt texture;
  TextureVertex()
  {
    position = new pt();
    texture = new pt();
  }
}
class ImageTexture
{

  //int size;
  TextureVertex [][] coordinates;


  //float image_sizex, image_sizey;

  PImage texture;

  ImageTexture() {
  };

  void Init()
  {
   // size = 51;
    texture = loadImage("data/image001.png");
    coordinates = new TextureVertex[ImageArraySize][ImageArraySize];
    pt startPoint = new pt(200, 200);
    float maxImageSize = 225;
    float delta = maxImageSize / (ImageArraySize - 1);
    float delta_texture = 1.0f / (ImageArraySize - 1);

    for (int i = 0; i < ImageArraySize; i++)
      for (int j = 0; j < ImageArraySize; j++)
      {
        coordinates[i][j] = new TextureVertex();
        coordinates[i][j].position = new pt(startPoint.x + j * delta, startPoint.y + i* delta);
        coordinates[i][j].texture = new pt(j * delta_texture, i* delta_texture);
      }
  }

  void drawImage()
  {
    
    // test
   /* for (int i = 0; i < ImageArraySize; i++)
      for (int j = 0; j < ImageArraySize; j++)
      {
        vertex(coordinates[i][j].position.x,coordinates[i][j].position.y);
      }
      */
    for (int i = 0; i < ImageArraySize - 1; i++)
    {
      noStroke();
      textureMode(NORMAL);
      beginShape();
      texture(texture);
        for (int j = 0; j < ImageArraySize - 1; j++)
        {
          vertex(coordinates[i][j].position.x, coordinates[i][j].position.y, coordinates[i][j].texture.x, coordinates[i][j].texture.y);
          vertex(coordinates[i + 1][j].position.x, coordinates[i + 1][j].position.y, coordinates[i + 1][j].texture.x, coordinates[i + 1][j].texture.y);
          vertex(coordinates[i + 1][j + 1].position.x, coordinates[i + 1][j + 1].position.y, coordinates[i + 1][j + 1].texture.x, coordinates[i + 1][j + 1].texture.y);
          vertex(coordinates[i][j + 1].position.x, coordinates[i][j + 1].position.y, coordinates[i][j + 1].texture.x, coordinates[i][j + 1].texture.y);
        }
          
     



      endShape();
    }
  }
}