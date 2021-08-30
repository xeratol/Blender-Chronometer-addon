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

def bridge_upper_lower_teeth(numVerts, startIdxUpper, startIdxLower):
    faces = []
    for i in range(numVerts - 1):
        face = (i + startIdxUpper + 1, i + startIdxUpper, i + startIdxLower, i + startIdxLower + 1)
        faces.append(face)
    face = (startIdxUpper, startIdxUpper + numVerts - 1, startIdxLower + numVerts - 1, startIdxLower)
    faces.append(face)
    return faces

def add_faces(numVertsTeeth, vertPerTooth,
        startIdxUpperTeeth, endIdxUpperTeeth,
        startIdxLowerTeeth, endIdxLowerTeeth):
        #startIdxUpperBase, endIdxUpperBase,
        #startIdxLowerBase, endIdxLowerBase):
    faces = []

    faces.extend(
        bridge_upper_lower_teeth( numVertsTeeth,
            startIdxUpperTeeth, startIdxLowerTeeth
        )
    )

    return faces

def create_teeth(vertPerTooth, numSegments, radius, addendum, z):
    verts = []
    for i in range(numSegments):
        angleRad = math.radians(i * 360.0 / numSegments)
        if (i % vertPerTooth == 0):
            radiusAdd = addendum
            vert = polar_coords(radius + radiusAdd, angleRad, z)
            verts.append(vert)
        radiusAdd = ( ( (i % vertPerTooth ) / (vertPerTooth) ) ** 4 ) * addendum
        vert = polar_coords(radius + radiusAdd, angleRad, z)
        verts.append(vert)
    return verts

def create_base(radius, numSegments, z):
    angleRad = math.radians( 360.0 / numSegments)
    verts = [polar_coords(radius, angleRad * i, z) for i in range(numSegments) ]
    return verts

def add_escape_wheel(self, context):
    verts = []
    
    numSegments = (self.vertPerTooth - 1) * self.numTeeth
    
    vertsUpperTeeth = create_teeth(self.vertPerTooth - 1, numSegments, self.radius, self.addendum, self.width / 2.0)
    vertsLowerTeeth = create_teeth(self.vertPerTooth - 1, numSegments, self.radius, self.addendum, -self.width / 2.0)
    vertsUpperTeethStartIdx = len(verts)
    verts.extend(vertsUpperTeeth)
    vertsLowerTeethStartIdx = len(verts)
    verts.extend(vertsLowerTeeth)

    vertsUpperBase = create_base(self.radius - self.base, numSegments, self.width / 2.0)
    vertsLowerBase = create_base(self.radius - self.base, numSegments, -self.width / 2.0)
    vertsUpperBaseStartIdx = len(verts)
    verts.extend(vertsUpperBase)
    vertsLowerBaseStartIdx = len(verts)
    verts.extend(vertsLowerBase)

    faces = add_faces(vertsLowerTeethStartIdx, self.vertPerTooth,
        vertsUpperTeethStartIdx, vertsLowerTeethStartIdx - 1,
        vertsLowerTeethStartIdx, vertsUpperBaseStartIdx - 1,
        vertsUpperBaseStartIdx, vertsLowerBaseStartIdx - 1,
        vertsLowerBaseStartIdx, len(verts) - 1)

    mesh = bpy.data.meshes.new(name="Chronometer")
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

    tanDist: FloatProperty(
        name="Tangential Distance",
        description="Tangential Distance of the Impulse Roller",
        min=0.0,
        soft_max=1000.0,
        unit='LENGTH',
        default=5.0
    )
    
    addendum: FloatProperty(
        name="Addendum",
        description="Addendum, extent of tooth above radius",
        min=0.0,
        soft_max=1000.0,
        unit='LENGTH',
        default=2.0
    )
    
    base: FloatProperty(
        name="Base",
        description="Base, extent of gear below radius",
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

        box = layout.box()
        box.prop(self, 'numTeeth')
        box.prop(self, 'vertPerTooth')

        box = layout.box()
        box.prop(self, 'radius')
        box.prop(self, 'tanDist')
        box.prop(self, 'base')
        box.prop(self, 'width')
        box.prop(self, 'addendum')

        box = layout.box()
        box.prop(self, 'internality')

        box = layout.box()
        box.prop(self, 'align', expand=True)
        box.prop(self, 'location', expand=True)
        box.prop(self, 'rotation', expand=True)
            
    def execute(self, context):

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