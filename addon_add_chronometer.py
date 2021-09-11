# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

bl_info = {
    "name": "Chronometer Escapement",
    "author": "Francis Joseph Serina",
    "version": (0, 0, 1),
    "blender": (2, 80, 0),
    "location": "View3D > Add > Mesh",
    "description": "Adds 5 pieces of the Chronometer Escapement - Recoil Type. The pieces are Escape Wheel, Impulse Roller, Spring of Detent, Discharge Pallet, and Locking Pallet.",
    "warning": "",
    "doc_url": "",
    "category": "Add Mesh",
}

import bpy
from bpy.types import Operator
from bpy.props import (
    FloatProperty,
    IntProperty,
)
from bpy_extras.object_utils import AddObjectHelper, object_data_add
import mathutils
import math

def polar_coords(radius, angleRad, z = 0):
    vert = (radius * math.cos(angleRad), radius * math.sin(angleRad), z)
    return vert

def flip_face(face):
    newFace = face[::-1]
    return newFace

def flip_faces(faces):
    newFaces = []
    for f in faces:
        newFaces.append(flip_face(f))
    return newFaces

def move_verts(verts, offset):
    for i, v in enumerate(verts):
        verts[i] = [v[x] + offset[x] for x in range(3)]

def rot_verts(verts, angleRad):
    cosTheta = math.cos( angleRad )
    sinTheta = math.sin( angleRad )

    for i, v in enumerate(verts):
        verts[i] = [v[0] * cosTheta - v[1] * sinTheta, v[0] * sinTheta + v[1] * cosTheta, v[2]]

def bridge_upper_lower_teeth(numVerts, startIdxUpper, startIdxLower):
    faces = []
    for i in range(numVerts - 1):
        face = (i + startIdxUpper + 1, i + startIdxUpper, i + startIdxLower, i + startIdxLower + 1)
        faces.append(face)
    face = (startIdxUpper, startIdxUpper + numVerts - 1, startIdxLower + numVerts - 1, startIdxLower)
    faces.append(face)
    return faces

def bridge_loops(numVerts, startIdxUpper, startIdxLower):
    faces = []
    for i in range(numVerts - 1):
        face = (i + startIdxUpper, i + startIdxUpper + 1, i + startIdxLower + 1, i + startIdxLower)
        faces.append(face)
    return faces

def bridge_teeth_base(vertPerTooth, baseStartIdx, baseEndIdx, teethStartIdx, teethEndIdx):
    faces = []
    j = teethStartIdx
    face = (baseStartIdx + 1, baseStartIdx, teethEndIdx, teethStartIdx, teethStartIdx + 1)
    faces.append(face)
    for i in range(baseStartIdx + 1, baseEndIdx):
        j += 1
        if ((j + 1) % vertPerTooth == 0):
            face = (i + 1, i, j, j + 1, j + 2)
            j += 1
        else:
            face = (i + 1, i, j, j + 1)
        faces.append(face)
    face = (baseStartIdx, baseEndIdx, teethEndIdx - 1, teethEndIdx)
    faces.append(face)
    return faces

def add_faces(numVertsTeeth, vertPerTooth,
        startIdxUpperTeeth, endIdxUpperTeeth,
        startIdxLowerTeeth, endIdxLowerTeeth,
        startIdxUpperBase, endIdxUpperBase,
        startIdxLowerBase, endIdxLowerBase):
    faces = []

    faces.extend(
        bridge_upper_lower_teeth( numVertsTeeth,
            startIdxUpperTeeth, startIdxLowerTeeth
        )
    )

    faces.extend(
        bridge_teeth_base( vertPerTooth,
            startIdxUpperBase, endIdxUpperBase,
            startIdxUpperTeeth, endIdxUpperTeeth
        )
    )

    faces.extend( flip_faces(
        bridge_teeth_base( vertPerTooth,
            startIdxLowerBase, endIdxLowerBase,
            startIdxLowerTeeth, endIdxLowerTeeth
        )
    ))

    return faces

def create_teeth(numTeeth, vertsPerTooth, radius, dedendum, teethShaperTheta, thetaPerTooth, z):
    teethShaperRadius = radius * math.sin( thetaPerTooth ) / math.sin( teethShaperTheta / 2 )
    teethShaperDist = radius + teethShaperRadius - dedendum
    toothSegmentAngleInc = teethShaperTheta * 0.5 / (vertsPerTooth - 1)

    verts = []
    for i in range(1, numTeeth + 1):
        toothAngle = thetaPerTooth * i
        toothVerts = []
        for j in range(vertsPerTooth):
            toothSegmentAngle = math.pi + toothAngle + j * toothSegmentAngleInc
            toothVerts.append( polar_coords(teethShaperRadius, toothSegmentAngle, z) )
        teethShaperCenter = polar_coords(teethShaperDist, toothAngle, 0)
        move_verts(toothVerts, teethShaperCenter)
        verts.extend( toothVerts[::-1] )
    return verts

def create_base(radius, numSegments, z):
    angleRad = 2 * math.pi / numSegments
    verts = [polar_coords(radius, angleRad * i, z) for i in range(numSegments) ]
    return verts

def create_arc(radius, numSegments, z, arc):
    angleRad = arc / (numSegments - 1)
    verts = [polar_coords(radius, angleRad * i, z) for i in range(numSegments) ]
    return verts

def add_escape_wheel(self, context):
    verts = []

    vertsUpperTeeth = create_teeth(self.numTeeth, self.vertPerTooth, self.radius, self.dedendum, self.teethShaperTheta, self.escWheelTheta, 0)
    vertsLowerTeeth = create_teeth(self.numTeeth, self.vertPerTooth, self.radius, self.dedendum, self.teethShaperTheta, self.escWheelTheta, -self.width)
    vertsUpperTeethStartIdx = len(verts)
    verts.extend(vertsUpperTeeth)
    vertsLowerTeethStartIdx = len(verts)
    verts.extend(vertsLowerTeeth)

    numSegments = (self.vertPerTooth - 1) * self.numTeeth
    base = self.radius - self.dedendum - self.escWheelBase

    vertsUpperBase = create_base(base, numSegments, 0)
    vertsLowerBase = create_base(base, numSegments, -self.width)
    vertsUpperBaseStartIdx = len(verts)
    verts.extend(vertsUpperBase)
    vertsLowerBaseStartIdx = len(verts)
    verts.extend(vertsLowerBase)

    rot_verts(verts, math.pi / 2)

    faces = add_faces(vertsLowerTeethStartIdx, self.vertPerTooth,
        vertsUpperTeethStartIdx, vertsLowerTeethStartIdx - 1,
        vertsLowerTeethStartIdx, vertsUpperBaseStartIdx - 1,
        vertsUpperBaseStartIdx, vertsLowerBaseStartIdx - 1,
        vertsLowerBaseStartIdx, len(verts) - 1)

    mesh = bpy.data.meshes.new(name="Escape Wheel")
    mesh.from_pydata(verts, [], faces)
    # useful for development when the mesh may be invalid.
    mesh.validate(verbose=True)
    obj = object_data_add(context, mesh, operator=self)

    vertGrp = obj.vertex_groups.new(name="Upper Teeth")
    vertGrp.add(list(range(vertsUpperTeethStartIdx, vertsLowerTeethStartIdx)), 1.0, 'ADD')
    vertGrp = obj.vertex_groups.new(name="Lower Teeth")
    vertGrp.add(list(range(vertsLowerTeethStartIdx, vertsUpperBaseStartIdx)), 1.0, 'ADD')
    vertGrp = obj.vertex_groups.new(name="Upper Base")
    vertGrp.add(list(range(vertsUpperBaseStartIdx, vertsLowerBaseStartIdx)), 1.0, 'ADD')
    vertGrp = obj.vertex_groups.new(name="Lower Base")
    vertGrp.add(list(range(vertsLowerBaseStartIdx, len(verts))), 1.0, 'ADD')

def add_impulse_roller(self, context):
    biggerArc = 2 * math.pi - self.impRollerTheta

    verts = []
    startIdxUpperOuter = len(verts)
    verts.extend( create_arc(self.impRollerRadius, self.impRollerVert, 0, biggerArc) )
    startIdxLowerOuter = len(verts)
    verts.extend( create_arc(self.impRollerRadius, self.impRollerVert, -self.width, biggerArc) )

    base = max(0, self.impRollerRadius - self.impRollerBase)
    startIdxUpperInner = len(verts)
    verts.extend( create_arc(base, self.impRollerVert, 0, biggerArc) )
    startIdxLowerInner = len(verts)
    verts.extend( create_arc(base, self.impRollerVert, -self.width, biggerArc) )

    rot_verts(verts, math.pi + math.atan( self.impRollerCenter[1] / self.impRollerCenter[0] ) - self.impRollerTheta / 2)

    faces = []
    faces.extend( bridge_loops(self.impRollerVert, startIdxLowerOuter, startIdxUpperOuter) )
    # faces.extend( bridge_loops(self.impRollerVert, startIdxLowerInner, startIdxUpperInner) )
    faces.extend( bridge_loops(self.impRollerVert, startIdxUpperOuter, startIdxUpperInner) )
    faces.extend( bridge_loops(self.impRollerVert, startIdxLowerInner, startIdxLowerOuter) )
    faces.append( [startIdxLowerOuter, startIdxUpperOuter, startIdxUpperInner, startIdxLowerInner] )
    faces.append( [len(verts) - 1, startIdxLowerInner - 1, startIdxLowerOuter - 1, startIdxUpperInner - 1] )

    mesh = bpy.data.meshes.new(name="Impulse Roller")
    mesh.from_pydata(verts, [], faces)
    # useful for development when the mesh may be invalid.
    mesh.validate(verbose=True)
    obj = object_data_add(context, mesh, operator=self)
    obj.location = self.impRollerCenter

    vertGrp = obj.vertex_groups.new(name="Upper Outer")
    vertGrp.add(list(range(startIdxUpperOuter, startIdxLowerOuter)), 1.0, 'ADD')

    vertGrp = obj.vertex_groups.new(name="Lower Outer")
    vertGrp.add(list(range(startIdxLowerOuter, startIdxUpperInner)), 1.0, 'ADD')

    vertGrp = obj.vertex_groups.new(name="Upper Inner")
    vertGrp.add(list(range(startIdxUpperInner, startIdxLowerInner)), 1.0, 'ADD')

    vertGrp = obj.vertex_groups.new(name="Lower Inner")
    vertGrp.add(list(range(startIdxLowerInner, len(verts))), 1.0, 'ADD')

def add_detent(self, context):
    verts = []

    verts.append( (0, 1, 0) )
    verts.append( (0, 0, 0) )
    verts.append( (0, 0, self.detentWidth) )
    verts.append( (0, 1, self.detentWidth) )

    verts.append( (self.detentLength, 1, 0) )
    verts.append( (self.detentLength, 0, 0) )
    verts.append( (self.detentLength, 0, self.detentWidth) )
    verts.append( (self.detentLength, 1, self.detentWidth) )

    faces = []

    faces.append( (0, 1, 2, 3) )
    faces.append( (0, 3, 7, 4) )
    faces.append( (3, 2, 6, 7) )
    faces.append( (2, 1, 5, 6) )
    faces.append( (1, 0, 4, 5) )
    faces.append( (4, 7, 6, 5) )

    mesh = bpy.data.meshes.new(name="Detent")
    mesh.from_pydata(verts, [], faces)
    # useful for development when the mesh may be invalid.
    mesh.validate(verbose=True)

    obj = object_data_add(context, mesh, operator=self)
    obj.location = (self.detentLeftEnd + self.impRollerCenter[0], self.impRollerCenter[1], 0)

def add_locking_pallet(self, context):
    verts = []

    verts.append( (0, 1, -self.width) )
    verts.append( (0, -self.lockingPalletDepth, -self.width) )
    verts.append( (0, -self.lockingPalletDepth, 0) )
    verts.append( (0, 1, 0) )
    
    verts.append( (1, 1, -self.width) )
    verts.append( (1, -self.lockingPalletDepth + 1, -self.width) )
    verts.append( (1, -self.lockingPalletDepth + 1, 0) )
    verts.append( (1, 1, 0) )

    move_verts(verts, (-self.detentLeftEnd - self.impRollerCenter[0], 0, 0))

    faces = []

    faces.append( (0, 1, 2, 3) )
    faces.append( (0, 3, 7, 4) )
    faces.append( (3, 2, 6, 7) )
    faces.append( (2, 1, 5, 6) )
    faces.append( (1, 0, 4, 5) )
    faces.append( (4, 5, 6, 7) )

    mesh = bpy.data.meshes.new(name="Locking Pallet")
    mesh.from_pydata(verts, [], faces)
    # useful for development when the mesh may be invalid.
    mesh.validate(verbose=True)

    obj = object_data_add(context, mesh, operator=self)
    obj.location = (self.detentLeftEnd + self.impRollerCenter[0], self.impRollerCenter[1], 0)


def add_discharge_pallet(self, context):
    verts = []

    verts.append( (self.dischargePalletTip[0], self.dischargePalletTip[1], self.detentWidth) )

    faces = []

    mesh = bpy.data.meshes.new(name="Discharge Pallet")
    mesh.from_pydata(verts, [], faces)
    # useful for development when the mesh may be invalid.
    mesh.validate(verbose=True)

    obj = object_data_add(context, mesh, operator=self)
    obj.location = self.impRollerCenter


class AddChronometer(Operator, AddObjectHelper):
    """Create a new Chronometer Escapement"""
    bl_idname = "mesh.primitive_chronometer"
    bl_label = "Add Chronometer Escapement"
    bl_options = {'REGISTER', 'UNDO', 'PRESET'}

    # Escape Wheel Properties

    numTeeth: IntProperty(
        name="Number of Teeth",
        description="Number of teeth of the Escape Wheel",
        default=12,
        min=9,
        soft_max=1000,
    )

    vertPerTooth: IntProperty(
        name="Vertices per Tooth",
        description="Number of Vertices per tooth, more for smoother",
        default=8,
        min=6,
        soft_max=100
    )
    
    radius: FloatProperty(
        name="Radius",
        description="Radius of the Escape Wheel",
        min=10.0,
        soft_max=1000.0,
        unit='LENGTH',
        default=10.0
    )

    dedendum: FloatProperty(
        name="Dedendum",
        description="Dedendum, extent of tooth below radius",
        min=0.0,
        soft_max=1000.0,
        unit='LENGTH',
        default=2.5
    )
    
    escWheelBase: FloatProperty(
        name="Base",
        description="Base, extent of escape wheel below radius",
        min=0.0,
        soft_max=1000.0,
        unit='LENGTH',
        default=2.0
    )
    
    escWheelTheta = 0.0
    
    teethShaperTheta: FloatProperty(
        name="Teeth Shaper Theta",
        description="Angle of Teeth Shaper overlapping with Escape Wheel (hidden)"
    )
    
    # Impulse Roller Properties

    restTooth: IntProperty(
        name="Resting Tooth",
        description="Tooth of the Escape Wheel where the Impulse Roller rests on",
        min=1,
        soft_max=100,
        default=1
    )
    
    impRollerVert: IntProperty(
        name="Vertices",
        description="Vertices of the Impulse Roller",
        default=32,
        min=6,
        soft_max=1000,
    )
    
    impRollerBase: FloatProperty(
        name="Base",
        description="Base, extent of impulse roller below radius",
        min=0.0,
        soft_max=1000.0,
        unit='LENGTH',
        default=2.0
    )

    impRollerCenter = [0.0, 0.0]
    impRollerRadius = 0.0
    impRollerTheta = 0.0

    # Detent Properties

    lockingPalletDepth: FloatProperty(
        name="Depth of the Locking Pallet",
        description="Depth of the Pallet locking the Escape Wheel",
        min=1.0,
        max=100.0,
        default=1.0,
        unit='LENGTH',
    )

    dischargeAngle: FloatProperty(
        name="Angle of Discharge Pallet",
        description="Angle of Discharge Pallet for Escape Wheel to release energy",
        min=0.0,
        max=math.pi / 2,
        default=math.pi / 4,
        step=1,
        precision=2,
        subtype='ANGLE',
    )

    dischargeRadius: FloatProperty(
        name="Discharge Pallet Radius Factor",
        description="Radius of Discharge Pallet relative to Radius of Impulse Roller",
        min=1.0,
        max=100.0,
        default=100.0,
        step=10,
        precision=1,
        subtype='PERCENTAGE',
    )

    detentWidth: FloatProperty(
        name="Width",
        description="Width, thickness of the Detent",
        min=1.0,
        soft_max=100.0,
        unit='LENGTH',
        default=5.0
    )

    dischargePalletTip = [0.0, 0.0]
    detentLeftEnd = 0.0
    detentLength = 0.0

    # Common Properties

    width: FloatProperty(
        name="Width",
        description="Width, thickness of Escape Wheel and Impulse Roller",
        min=0.0,
        soft_max=1000.0,
        unit='LENGTH',
        default=5.0
    )

    def draw(self, context):
        layout = self.layout

        layout.label(text="Escape Wheel")
        box = layout.box()
        box.prop(self, 'numTeeth')
        box.prop(self, 'vertPerTooth')
        box.prop(self, 'radius')
        box.prop(self, 'escWheelBase')
        box.prop(self, 'dedendum')
        box.prop(self, 'width') # shared with Impulse Roller

        layout.label(text="Impulse Roller")
        box = layout.box()
        box.prop(self, 'restTooth')
        box.prop(self, 'impRollerVert')
        box.prop(self, 'impRollerBase')
        box.prop(self, 'width') # shared with Escape Wheel

        layout.label(text="Detent")
        box = layout.box()
        box.prop(self, 'lockingPalletDepth')
        box.prop(self, 'dischargeAngle')
        box.prop(self, 'dischargeRadius')
        box.prop(self, 'detentWidth')

        box = layout.box()
        box.prop(self, 'align', expand=True)
        box.prop(self, 'location', expand=True)
        box.prop(self, 'rotation', expand=True)


    def validate(self):
        self.escWheelTheta = 2 * math.pi / self.numTeeth

        # Escape Wheel Computations
        maxRestTooth = self.numTeeth // 4 - 1
        self.restTooth = min(maxRestTooth, self.restTooth)

        x = self.radius * ( 1 - math.cos( self.escWheelTheta ) ) - self.dedendum
        x /= self.radius * math.sin( self.escWheelTheta )

        if (x > -0.1):
            self.report({'WARNING'}, 'Invalid Escape Wheel dimensions. Dedendum automatically adjusted.')
            x = -0.1
            minDedendum = self.radius * ( 1 - math.cos(self.escWheelTheta) + math.sin(self.escWheelTheta) * 0.1 )
            self.dedendum = minDedendum

        self.teethShaperTheta = -4 * math.atan( x )

        if (self.teethShaperTheta > math.pi - 2 * self.escWheelTheta):
            self.report({'WARNING'}, 'Invalid Teeth dimensions. Dedendum automatically adjusted.')
            self.teethShaperTheta = math.pi - 2 * self.escWheelTheta
            maxDedendum = self.radius * ( 1 - math.cos(self.escWheelTheta) + math.sin(self.escWheelTheta) * math.tan(self.teethShaperTheta / 4) )
            self.dedendum = maxDedendum

        if (self.radius - self.dedendum <= self.escWheelBase):
            self.report({'WARNING'}, 'Invalid Escape Wheel Dimensions. Base automatically adjusted.')
            self.escWheelBase = self.radius - self.dedendum # maximum escape wheel base

        # Impulse Roller Computations
        tanDist = self.radius / math.tan( ( math.pi - self.escWheelTheta * (2 * self.restTooth + 1) ) / 2 )
        self.impRollerCenter = [tanDist, self.radius, 0]
        distBetween = math.sqrt( self.impRollerCenter[0] ** 2 + self.impRollerCenter[1] ** 2 )
        self.impRollerTheta = 2 * math.atan( ( self.impRollerCenter[1] * math.sin( self.escWheelTheta / 2 ) ) / ( distBetween - self.impRollerCenter[1] * math.cos( self.escWheelTheta / 2 )) )
        self.impRollerRadius = ( self.impRollerCenter[1] * math.sin( self.escWheelTheta / 2 ) ) / math.sin( self.impRollerTheta / 2 )

        # Detent Computations
        i = [ -self.impRollerCenter[0], self.lockingPalletDepth ]
        r = self.impRollerRadius * self.dischargeRadius / 100
        minLeftEnd = -self.impRollerCenter[0] - 2 * self.radius
        maxLeftEnd = -self.impRollerCenter[0] - 0.1 * self.radius
        m = math.tan( math.asin( r / abs( minLeftEnd ) ) )
        maxIY = min( m * ( i[0] - minLeftEnd ), self.dedendum)
        if (i[1] > maxIY):
            self.lockingPalletDepth = i[1] = maxIY
        #self.report({'INFO'}, "{0}".format(i[1]))

        if ( (self.impRollerCenter[0] - i[0]) ** 2 + (self.impRollerCenter[1] - i[1]) ** 2 >= (self.impRollerRadius * 1.2) ** 2 ):
            r = self.impRollerRadius * self.dischargeRadius / 100
            self.dischargePalletTip = [
                r * math.cos( math.pi - self.dischargeAngle ),
                r * math.sin( self.dischargeAngle )
            ]

            minDetentAngle = math.atan2( i[1], (i[0] - minLeftEnd) )
            maxDetentAngle = math.atan2( i[1], (i[0] - maxLeftEnd) )
            #self.report({'INFO'}, '{0} - {1}'.format( math.degrees(minDetentAngle), math.degrees(maxDetentAngle)))
            detentAngle = math.atan2( self.dischargePalletTip[1] - i[1], self.dischargePalletTip[0] - i[0] )
            #self.report({'INFO'}, "{0}".format(math.degrees(detentAngle)))

            if (detentAngle > maxDetentAngle):
                detentAngle = maxDetentAngle
            elif (detentAngle < minDetentAngle):
                detentAngle = minDetentAngle
            m = math.tan( detentAngle )
            #self.report({'INFO'}, "{0}".format(math.degrees(detentAngle)))
            self.detentLeftEnd = i[0] - i[1] / m

            minDetentLength = abs( self.detentLeftEnd ) * math.cos( detentAngle ) - math.sqrt( self.impRollerRadius ** 2 - ( self.detentLeftEnd * math.sin( detentAngle ) ) ** 2 )
            minDischargeAngle = math.asin( minDetentLength * math.sin( detentAngle ) / self.impRollerRadius )
            maxDischargeAngle = ( math.pi * 0.5 ) - detentAngle
            #self.report({'INFO'}, '{0} - {1}'.format( math.degrees(minDischargeAngle), math.degrees(maxDischargeAngle)))

            if (self.dischargeAngle > maxDischargeAngle):
                self.dischargeAngle = maxDischargeAngle
            elif ( self.dischargeAngle < minDischargeAngle ):
                self.dischargeAngle = minDischargeAngle

            r = math.sin( detentAngle ) * abs( self.detentLeftEnd ) / math.sin( math.pi - detentAngle - self.dischargeAngle )
            self.dischargeRadius = max( 1, min( 100 * r / self.impRollerRadius, 100 ) )

            r = self.impRollerRadius * self.dischargeRadius / 100
            self.dischargePalletTip = [
                r * math.cos( math.pi - self.dischargeAngle ),
                r * math.sin( self.dischargeAngle )
            ]

            self.detentLength = math.sqrt( ( self.detentLeftEnd - self.dischargePalletTip[0] ) ** 2 + ( self.dischargePalletTip[1] ) ** 2 )

        else:
            self.report({'ERROR'}, 'Not enough space between Impulse Roller and Locking Pallet')

        return True

            
    def execute(self, context):

        if (not self.validate()):
            return
            
        add_impulse_roller(self, context)
        add_escape_wheel(self, context)
        add_detent(self, context)
        add_locking_pallet(self, context)
        add_discharge_pallet(self, context)

        return {'FINISHED'}


# Registration

def add_object_button(self, context):
    self.layout.operator(
        AddChronometer.bl_idname,
        text="Chronometer Escapement",
        icon='PIVOT_INDIVIDUAL')


def register():
    bpy.utils.register_class(AddChronometer)
    bpy.types.VIEW3D_MT_mesh_add.append(add_object_button)


def unregister():
    bpy.utils.unregister_class(AddChronometer)
    bpy.types.VIEW3D_MT_mesh_add.remove(add_object_button)


if __name__ == "__main__":
    register()
