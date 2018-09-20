import cv2

# Windows dependencies
# - Python 2.7.6: http://www.python.org/download/
# - OpenCV: http://opencv.org/
# - Numpy -- get numpy from here because the official builds don't support x64:
#   http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy

# Mac Dependencies
# - brew install python
# - pip install numpy
# - brew tap homebrew/science
# - brew install opencv

cap = cv2.VideoCapture(0)
cap.set(3, 640)
cap.set(4, 480)
cap.set(15, 0.1)

while(True):
    ret, frame = cap.read()
    rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2BGRA)


    cv2.imshow('frame', rgb)
    if cv2.waitKey(1) & 0xFF == ord('q'):
        newframe = cv2.resize(frame,(28,28))
        out = cv2.imwrite('capture.jpg', frame)
        newout = cv2.imwrite('newcapture.jpg', newframe)

        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        ret, gb = cv2.threshold(gray,128,255,cv2.THRESH_BINARY)
        gb = cv2.bitwise_not(gb)
        newgb = cv2.resize(gb,(28,28))
        newgbout = cv2.imwrite('newgb.jpg', newgb)
        break

cap.release()
cv2.destroyAllWindows()