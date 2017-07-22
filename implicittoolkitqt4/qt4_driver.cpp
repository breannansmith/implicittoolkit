#include <QApplication>
#include <QDesktopWidget>

#include "Window.h"

int main( int argc, char **argv )
{
  QApplication app(argc, argv);
  Window window;
  window.resize(window.sizeHint());
  window.setWindowTitle("SDF Computer");
  window.show();
  window.raise();
  return app.exec();
}
