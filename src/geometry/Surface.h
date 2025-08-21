#pragma once

class Surface
{
public:
  Vector3d refpt, normdir;
  virtual void display(const std::vector<Vector3d>& vertices);
  virtual void reverse(void);
  virtual int operator==(const Surface& other);
  virtual int pointMember(std::vector<Vector3d>& vertices, Vector3d pt);
};

class CylinderSurface : public Surface
{
public:
  CylinderSurface(Vector3d center, Vector3d normdir, double r);
  void display(const std::vector<Vector3d>& vertices);
  void reverse(void);
  int operator==(const CylinderSurface& other);
int operator==(const Surface& other) override {
    // Optionally, dynamic_cast and call the CylinderSurface overload
    if (auto arc = dynamic_cast<const CylinderSurface*>(&other)) {
        return (*this == *arc);
    }
    return 0;
}

virtual int pointMember(std::vector<Vector3d>& vertices, Vector3d pt);

  double r;
};
