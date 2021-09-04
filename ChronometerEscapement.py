# References:
# https://www.youtube.com/watch?v=nJ_k4VbAeGY
# https://en.wikipedia.org/wiki/Escapement#/media/File:Diagram_of_the_chronometer_detent_escapement_Britten's_Clock_and_Watchmaker's_Handbook_9th_Edition_1896.jpg
# https://www.gear2motion.com/lever

from PIL import Image, ImageDraw
import math

color = (255, 255, 255)
imageDimensions = (2000, 2000)
imageCenter = (imageDimensions[0]//2, imageDimensions[1]//2)

class EscapeWheel:
    center = [0, 0]
    radius = 500
    n = 6 # must be >= 4

    def getTheta(self): # radians
        return 2 * math.pi / self.n

    def getAbsCenter(self):
        return (imageCenter[0] + self.center[0], imageCenter[1] + self.center[1])

    def getLowerPoint(self):
        absCenter = self.getAbsCenter()
        return (absCenter[0] - self.radius, absCenter[1] - self.radius)

    def getHigherPoint(self):
        absCenter = self.getAbsCenter()
        return (absCenter[0] + self.radius, absCenter[1] + self.radius)

    def drawCircle(self, painter):
        painter.ellipse( [self.getLowerPoint(), self.getHigherPoint()], fill=None, outline=color )

    def getVertices(self, angleRadOffset = 0):
        absCenter = self.getAbsCenter()
        angleInc = 2 * math.pi / self.n
        verts = []
        for i in range(self.n):
            x = absCenter[0] + self.radius * math.cos(i * angleInc + angleRadOffset)
            y = absCenter[0] + self.radius * math.sin(i * angleInc + angleRadOffset)
            verts.append((x, y))
        return verts

    def drawPolygon(self, painter : ImageDraw, rotDeg = 0):
        painter.polygon( self.getVertices(angleRadOffset=math.radians(rotDeg)), outline=(255, 64, 64) )

class ImpulseRoller:
    center = [0, 0]
    radius = 500
    theta = 0 # radians

    def getAbsCenter(self):
        return (imageCenter[0] + self.center[0], imageCenter[1] + self.center[1])

    def getLowerPoint(self):
        absCenter = self.getAbsCenter()
        return (absCenter[0] - self.radius, absCenter[1] - self.radius)

    def getHigherPoint(self):
        absCenter = self.getAbsCenter()
        return (absCenter[0] + self.radius, absCenter[1] + self.radius)

    def drawCircle(self, painter : ImageDraw):
        painter.ellipse( [self.getLowerPoint(), self.getHigherPoint()], fill=None, outline=color )

    def draw(self, painter : ImageDraw, escapeWheel : EscapeWheel):
        painter.arc( [ escapeWheel.getLowerPoint(), escapeWheel.getHigherPoint() ],
            start=math.degrees( math.atan( self.center[1] / self.center[0] ) - ( escapeWheel.getTheta() / 2 ) ),
            end=math.degrees( math.atan( self.center[1] / self.center[0] ) + ( escapeWheel.getTheta() / 2 ) ),
            fill=(128, 128, 255) )

        painter.arc( [ self.getLowerPoint(), self.getHigherPoint() ],
            start=-90 - math.degrees( math.atan( self.center[0] / self.center[1] ) - ( self.theta / 2 ) ),
            end=-90 - math.degrees( math.atan( self.center[0] / self.center[1] ) + ( self.theta / 2 ) ),
            fill=(128, 128, 255) )

canvas = Image.new('RGB', imageDimensions)
painter = ImageDraw.Draw(canvas)

escapeWheel = EscapeWheel()
impulseRoller = ImpulseRoller()

# Input
escapeWheel.radius = 400
escapeWheel.n = 9
impulseRoller.center = (500, escapeWheel.radius)
dedendum = 200

distBetween = math.sqrt( (escapeWheel.center[0] - impulseRoller.center[0]) ** 2 + (escapeWheel.center[1] - impulseRoller.center[1]) ** 2 )

impulseRoller.theta = 2 * math.atan( ( escapeWheel.radius * math.sin( escapeWheel.getTheta() / 2 ) ) / ( distBetween - escapeWheel.radius * math.cos( escapeWheel.getTheta() / 2 )) )
impulseRoller.radius = ( escapeWheel.radius * math.sin( escapeWheel.getTheta() / 2 ) ) / math.sin( impulseRoller.theta / 2 )

if (impulseRoller.center[0] <= impulseRoller.radius or
    impulseRoller.center[1] <= impulseRoller.radius):
    print("Warning: Not enough space between Locking Pallet and Impulse Roller.")

# escapeWheel.drawCircle(painter)
escapeWheelAngleDegOffset = math.degrees( math.atan( impulseRoller.center[1] / impulseRoller.center[0] ) - ( escapeWheel.getTheta() / 2 ) )
escapeWheel.drawPolygon(painter, escapeWheelAngleDegOffset)

# impulseRoller.drawCircle(painter)
impulseRoller.draw(painter, escapeWheel)

# teeth shaper
x = escapeWheel.radius * ( 1 - math.cos( escapeWheel.getTheta() ) ) - dedendum
x /= escapeWheel.radius * math.sin( escapeWheel.getTheta() )
teethShaperTheta = -4 * math.atan( x )
# print(math.degrees(teethShaperTheta))
teethShaperRadius = escapeWheel.radius * math.sin( escapeWheel.getTheta() ) / math.sin( teethShaperTheta / 2 )
# print(teethShaperRadius)

if (teethShaperRadius < 0 or teethShaperTheta < 0):
    print("Warning: Cannot create Escape Wheel Teeth.")

distTeethToEscape = escapeWheel.radius + teethShaperRadius - dedendum
# print(distTeethToEscape)

for i in range(escapeWheel.n):
    angleRad = math.radians( escapeWheelAngleDegOffset ) + i * escapeWheel.getTheta()
    teethShaperCenter = [distTeethToEscape * math.cos(angleRad), distTeethToEscape * math.sin(angleRad)]
    # painter.line([escapeWheel.getAbsCenter(),
        # (escapeWheel.getAbsCenter()[0] + teethShaperCenter[0],
        # escapeWheel.getAbsCenter()[1] + teethShaperCenter[1])],
        # fill=(255,255,255))
    lowerPt = (escapeWheel.getAbsCenter()[0] + teethShaperCenter[0] - teethShaperRadius,
        escapeWheel.getAbsCenter()[1] + teethShaperCenter[1] - teethShaperRadius)
    higherPt = (escapeWheel.getAbsCenter()[0] + teethShaperCenter[0] + teethShaperRadius,
        escapeWheel.getAbsCenter()[1] + teethShaperCenter[1] + teethShaperRadius)
    if (teethShaperRadius < 0):
        boundingBox = [higherPt, lowerPt]
    else:
        boundingBox = [lowerPt, higherPt]
    painter.ellipse(boundingBox, fill=None, outline=(255,255,255))

# detent
painter.line([(imageCenter[0], (imageCenter[1] + escapeWheel.radius)),
    impulseRoller.getAbsCenter()],
    fill=(255, 255, 0))

canvas.show()