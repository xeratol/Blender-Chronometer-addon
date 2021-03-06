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
    "version": (0, 0, 2),
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
    verts = [
        (0, 1, 0),
        (0, 0, 0),
        (0, 0, self.detentWidth),
        (0, 1, self.detentWidth),
        (self.detentLength, 1, 0),
        (self.detentLength, 0, 0),
        (self.detentLength, 0, self.detentWidth),
        (self.detentLength, 1, self.detentWidth),
    ]

    faces = [
        (0, 1, 2, 3),
        (0, 3, 7, 4),
        (3, 2, 6, 7),
        (2, 1, 5, 6),
        (1, 0, 4, 5),
        (4, 7, 6, 5),
    ]

    mesh = bpy.data.meshes.new(name="Detent")
    mesh.from_pydata(verts, [], faces)
    # useful for development when the mesh may be invalid.
    mesh.validate(verbose=True)

    obj = object_data_add(context, mesh, operator=self)
    obj.location = (-self.detentBladeLength, self.impRollerCenter[1], 0)

def add_locking_pallet(self, context):
    verts = [
        (0, 1, -self.width),
        (0, -self.lockingPalletDepth, -self.width),
        (0, -self.lockingPalletDepth, 0),
        (0, 1, 0),
        (1, 1, -self.width),
        (1, -self.lockingPalletDepth + 1, -self.width),
        (1, -self.lockingPalletDepth + 1, 0),
        (1, 1, 0),
    ]
    move_verts(verts, (self.detentBladeLength, 0, 0))

    faces = [
        (0, 1, 2, 3),
        (0, 3, 7, 4),
        (3, 2, 6, 7),
        (2, 1, 5, 6),
        (1, 0, 4, 5),
        (4, 5, 6, 7),
    ]

    mesh = bpy.data.meshes.new(name="Locking Pallet")
    mesh.from_pydata(verts, [], faces)
    # useful for development when the mesh may be invalid.
    mesh.validate(verbose=True)

    obj = object_data_add(context, mesh, operator=self)
    obj.location = (-self.detentBladeLength, self.impRollerCenter[1], 0)

def add_discharge_pallet(self, context):
    adds = [[0,0]]
    theta = math.pi / 2 - self.detentReleaseAngle
    adds.append( [math.cos(theta), math.sin(-theta)] )
    theta = math.pi / 2 - self.dischargeAngle - self.detentReleaseAngle
    length = math.cos(theta)
    adds.append( [length * math.cos(-self.dischargeAngle), length * math.sin(-self.dischargeAngle)] )

    verts = [
        (self.dischargePalletTip[0] + adds[0][0], self.dischargePalletTip[1] + adds[0][1], 0),
        (self.dischargePalletTip[0] + adds[1][0], self.dischargePalletTip[1] + adds[1][1], 0),
        (self.dischargePalletTip[0] + adds[2][0], self.dischargePalletTip[1] + adds[2][1], 0),
        (self.dischargePalletTip[0] + adds[0][0], self.dischargePalletTip[1] + adds[0][1], self.detentWidth),
        (self.dischargePalletTip[0] + adds[1][0], self.dischargePalletTip[1] + adds[1][1], self.detentWidth),
        (self.dischargePalletTip[0] + adds[2][0], self.dischargePalletTip[1] + adds[2][1], self.detentWidth),
    ]

    faces = [
        (0, 2, 1),
        (0, 3, 5, 2),
        (1, 2, 5, 4),
        (0, 1, 4, 3),
        (3, 4, 5),
    ]

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

    impRollerCenter = [0.0, 0.0]
    impRollerRadius = 0.0
    impRollerTheta = 0.0
    dischargePalletTip = [0.0, 0.0] # relative to Impulse Roller Center
    dischargeRadiusFactor = 1.0     # relative to Impulse Roller Radius

    # Detent Properties

    lockingPalletDepth: FloatProperty(
        name="Depth of the Locking Pallet",
        description="Depth of the Pallet locking the Escape Wheel",
        min=1.0,
        max=100.0,
        default=1.0,
        unit='LENGTH',
    )

    detentBladeLength: FloatProperty(
        name="Blade of Detent Length",
        description="Length of the Detent, to the left of the Locking Pallet",
        min=1.0,
        max=1000.0,
        default=10.0,
        unit='LENGTH',
    )

    detentWidth: FloatProperty(
        name="Width",
        description="Width, thickness of the Detent",
        min=1.0,
        soft_max=100.0,
        unit='LENGTH',
        default=5.0
    )

    detentLength = 0.0
    detentReleaseAngle = 0.0

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
        box.prop(self, 'dischargeAngle')
        box.prop(self, 'width') # shared with Escape Wheel

        layout.label(text="Detent")
        box = layout.box()
        box.prop(self, 'lockingPalletDepth')
        box.prop(self, 'detentBladeLength')
        box.prop(self, 'detentWidth')

        box = layout.box()
        box.prop(self, 'align', expand=True)
        box.prop(self, 'location', expand=True)
        box.prop(self, 'rotation', expand=True)


    def validate_escape_wheel(self):
        self.escWheelTheta = 2 * math.pi / self.numTeeth

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

        return True


    def validate_impulse_roller(self):
        tanDist = self.radius / math.tan( ( math.pi - self.escWheelTheta * (2 * self.restTooth + 1) ) / 2 )
        self.impRollerCenter = [tanDist, self.radius, 0]
        distBetween = math.sqrt( self.impRollerCenter[0] ** 2 + self.impRollerCenter[1] ** 2 )
        self.impRollerTheta = 2 * math.atan( ( self.impRollerCenter[1] * math.sin( self.escWheelTheta / 2 ) ) / ( distBetween - self.impRollerCenter[1] * math.cos( self.escWheelTheta / 2 )) )
        self.impRollerRadius = ( self.impRollerCenter[1] * math.sin( self.escWheelTheta / 2 ) ) / math.sin( self.impRollerTheta / 2 )

        chordLength = 2 * self.radius * math.sin( self.escWheelTheta / 2 )
        if (self.impRollerRadius * 2 < chordLength):
            self.report({'ERROR'}, 'Not enough space between Impulse Roller and Locking Pallet') # untested
            return False

        return True


    def validate_detent(self):
        # relative to self.impRollerCenter
        self.lockingPalletDepth = min( self.dedendum, self.lockingPalletDepth )

        releaseAngleMin = math.radians(1.0)
        #releaseAngleMax = math.radians(45.0)
        self.detentBladeLength = min( max( self.lockingPalletDepth, self.detentBladeLength ), self.lockingPalletDepth / math.tan( releaseAngleMin ) )

        self.detentLength = self.impRollerCenter[0] + self.detentBladeLength
        self.detentReleaseAngle = math.atan( self.lockingPalletDepth / self.detentBladeLength )

        dischargeAngleMin = math.pi - math.atan2( self.lockingPalletDepth, -self.impRollerCenter[0] )
        dischargeAngleMax = math.pi / 2 - self.detentReleaseAngle
        self.dischargeAngle = min( dischargeAngleMax, max(dischargeAngleMin, self.dischargeAngle) )

        dischargePalletAngle = math.pi - self.detentReleaseAngle - self.dischargeAngle
        dischargeRadius = self.detentLength * math.sin(self.detentReleaseAngle) / math.sin(dischargePalletAngle)
        self.dischargeRadiusFactor = dischargeRadius / self.impRollerRadius
        self.dischargePalletTip = [ dischargeRadius * math.cos( math.pi - self.dischargeAngle ), dischargeRadius * math.sin( self.dischargeAngle ) ]

        self.detentLength = math.sqrt( (self.dischargePalletTip[1]) ** 2 + (self.detentLength - abs(self.dischargePalletTip[0])) ** 2)

        return True


    def execute(self, context):

        if (not self.validate_escape_wheel() or not self.validate_impulse_roller() or not self.validate_detent()):
            return {'FINISHED'}
            
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
