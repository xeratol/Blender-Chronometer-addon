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
    "description": "Adds a new Chronometer Escapement",
    "warning": "",
    "doc_url": "",
    "category": "Add Mesh",
}

import bpy
from bpy.types import Operator
from bpy.props import (
    FloatProperty,
    IntProperty,
    EnumProperty
)
from bpy_extras.object_utils import AddObjectHelper, object_data_add
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
    for i in range(baseStartIdx, baseEndIdx):
        if (j % vertPerTooth == 0):
            face = (i + 1, i, j, j + 1, j + 2)
            j += 1
        else:
            face = (i + 1, i, j, j + 1)
        j += 1
        faces.append(face)
    face = (baseStartIdx, baseEndIdx, teethEndIdx - 1, teethEndIdx, teethStartIdx)
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

def create_teeth(vertPerTooth, numSegments, radius, dedendum, z):
    verts = []
    for i in range(numSegments):
        angleRad = math.radians(i * 360.0 / numSegments)
        if (i % vertPerTooth == 0):
            radiusAdd = 0 #dedendum
            vert = polar_coords(radius + radiusAdd, angleRad, z)
            verts.append(vert)
        radiusAdd = ( ( 1.0 - ( i % vertPerTooth ) / vertPerTooth ) ** 4 ) * dedendum
        vert = polar_coords(radius + radiusAdd, angleRad, z)
        verts.append(vert)
    return verts

def create_base(radius, numSegments, z):
    angleRad = math.radians( 360.0 / numSegments)
    verts = [polar_coords(radius, angleRad * i, z) for i in range(numSegments) ]
    return verts

def create_arc(radius, numSegments, z, arc):
    angleRad = arc / (numSegments - 1)
    verts = [polar_coords(radius, angleRad * i, z) for i in range(numSegments) ]
    return verts

def add_escape_wheel(self, context):
    numSegments = (self.vertPerTooth - 1) * self.numTeeth
    
    verts = []
    vertsUpperTeeth = create_teeth(self.vertPerTooth - 1, numSegments, self.radius - self.dedendum, self.dedendum, self.width / 2.0)
    vertsLowerTeeth = create_teeth(self.vertPerTooth - 1, numSegments, self.radius - self.dedendum, self.dedendum, -self.width / 2.0)
    vertsUpperTeethStartIdx = len(verts)
    verts.extend(vertsUpperTeeth)
    vertsLowerTeethStartIdx = len(verts)
    verts.extend(vertsLowerTeeth)

    base = max(0, self.radius - self.dedendum - self.escWheelBase)
    vertsUpperBase = create_base(base, numSegments, self.width / 2.0)
    vertsLowerBase = create_base(base, numSegments, -self.width / 2.0)
    vertsUpperBaseStartIdx = len(verts)
    verts.extend(vertsUpperBase)
    vertsLowerBaseStartIdx = len(verts)
    verts.extend(vertsLowerBase)

    impRollerCenter = [self.tanDist, self.radius, 0]
    escWheelTheta = 2 * math.pi / self.numTeeth
    angleRadOffset = math.atan( impRollerCenter[1] / impRollerCenter[0] ) - ( escWheelTheta / 2 )
    rot_verts(verts, angleRadOffset)

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
    center = [self.tanDist, self.radius, 0]
    distBetween = math.sqrt( center[0] ** 2 + center[1] ** 2 )
    escapeWheelTheta = 2 * math.pi / self.numTeeth
    impulseRollerTheta = 2 * math.atan( ( self.radius * math.sin( escapeWheelTheta / 2 ) ) / ( distBetween - self.radius * math.cos( escapeWheelTheta / 2 )) )
    impulseRollerRadius = ( self.radius * math.sin( escapeWheelTheta / 2 ) ) / math.sin( impulseRollerTheta / 2 )

    if (center[0] <= impulseRollerRadius or center[1] <= impulseRollerRadius):
        self.report({'WARNING'}, 'Not enough space between Locking Pallet and Impulse Roller.')

    biggerArc = 2 * math.pi - impulseRollerTheta

    verts = []
    startIdxUpperOuter = len(verts)
    verts.extend( create_arc(impulseRollerRadius, self.impRollerVert, 0, biggerArc) )
    startIdxLowerOuter = len(verts)
    verts.extend( create_arc(impulseRollerRadius, self.impRollerVert, -self.width / 2.0, biggerArc) )

    base = max(0, impulseRollerRadius - self.impRollerBase)
    startIdxUpperInner = len(verts)
    verts.extend( create_arc(base, self.impRollerVert, 0, biggerArc) )
    startIdxLowerInner = len(verts)
    verts.extend( create_arc(base, self.impRollerVert, -self.width / 2.0, biggerArc) )

    rot_verts(verts, math.pi + math.atan( center[1] / center[0] ) + impulseRollerTheta / 2)
    move_verts(verts, center)

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

    vertGrp = obj.vertex_groups.new(name="Upper Outer")
    vertGrp.add(list(range(startIdxUpperOuter, startIdxLowerOuter)), 1.0, 'ADD')

    vertGrp = obj.vertex_groups.new(name="Lower Outer")
    vertGrp.add(list(range(startIdxLowerOuter, startIdxUpperInner)), 1.0, 'ADD')

    vertGrp = obj.vertex_groups.new(name="Upper Inner")
    vertGrp.add(list(range(startIdxUpperInner, startIdxLowerInner)), 1.0, 'ADD')

    vertGrp = obj.vertex_groups.new(name="Lower Inner")
    vertGrp.add(list(range(startIdxLowerInner, len(verts))), 1.0, 'ADD')


class AddChronometer(Operator, AddObjectHelper):
    """Create a new Chronometer Escapement"""
    bl_idname = "mesh.primitive_chronometer"
    bl_label = "Add Chronometer Escapement"
    bl_options = {'REGISTER', 'UNDO', 'PRESET'}

    numTeeth: IntProperty(
        name="Number of Teeth",
        description="Number of teeth of the Escape Wheel",
        default=12,
        min=6,
        soft_max=1000,
    )
    
    vertPerTooth: IntProperty(
        name="Vertices per Tooth",
        description="Number of Vertices per tooth, more for smoother",
        default=4,
        min=2,
        soft_max=100
    )
    
    radius: FloatProperty(
        name="Radius",
        description="Radius of the Escape Wheel",
        min=0.0,
        soft_max=1000.0,
        unit='LENGTH',
        default=10.0
    )

    impRollerVert: IntProperty(
        name="Vertices",
        description="Vertices of the Impulse Roller",
        default=32,
        min=6,
        soft_max=1000,
    )
    
    tanDist: FloatProperty(
        name="Tangential Distance",
        description="Tangential Distance of the Impulse Roller to the Escape Wheel",
        min=0.1,
        soft_max=1000.0,
        unit='LENGTH',
        default=5.0
    )
    
    dedendum: FloatProperty(
        name="Dedendum",
        description="Dedendum, extent of tooth below radius",
        min=0.0,
        soft_max=1000.0,
        unit='LENGTH',
        default=2.0
    )
    
    escWheelBase: FloatProperty(
        name="Base",
        description="Base, extent of escape wheel below radius",
        min=0.0,
        soft_max=1000.0,
        unit='LENGTH',
        default=2.0
    )
    
    impRollerBase: FloatProperty(
        name="Base",
        description="Base, extent of impulse roller below radius",
        min=0.0,
        soft_max=1000.0,
        unit='LENGTH',
        default=2.0
    )

    width: FloatProperty(
        name="Width",
        description="Width, thickness of Escape Wheel",
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

        layout.label(text="Impulse Roller")
        box = layout.box()
        box.prop(self, 'tanDist')
        box.prop(self, 'impRollerVert')
        box.prop(self, 'impRollerBase')

        layout.label(text="Detent")

        layout.label(text="Common")
        box = layout.box()
        box.prop(self, 'width')

        box = layout.box()
        box.prop(self, 'align', expand=True)
        box.prop(self, 'location', expand=True)
        box.prop(self, 'rotation', expand=True)
            
    def execute(self, context):

        add_impulse_roller(self, context)
        add_escape_wheel(self, context)

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
