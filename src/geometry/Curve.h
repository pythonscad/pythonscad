#pragma once

class Curve
{
public:
  int start, end;
  virtual void display(const std::vector<Vector3d>& vertices);
  virtual void reverse(void);
  virtual int operator==(const Curve& other);
  virtual int pointMember(std::vector<Vector3d>& vertices, Vector3d pt);
};

class ArcCurve : public Curve
{
public:
  ArcCurve(Vector3d center, Vector3d normdir, double r);
  void display(const std::vector<Vector3d>& vertices);
  void reverse(void);
  int operator==(const ArcCurve& other);

  int operator==(const Curve& other) override {
    // Optionally, dynamic_cast and call the ArcCurve overload
    if (auto arc = dynamic_cast<const ArcCurve*>(&other)) {
        return (*this == *arc);
    }
    return 0;
}
  double calcAngle(Vector3d refdir, Vector3d dir, Vector3d normdir);
  virtual int pointMember(std::vector<Vector3d>& vertices, Vector3d pt);

  double r;
  Vector3d center, normdir;
};
