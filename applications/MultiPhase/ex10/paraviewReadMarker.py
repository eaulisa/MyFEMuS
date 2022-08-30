# trace generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
marker = CSVReader(registrationName='marker*', FileName=['/home/eaulisa/FEMuS/femusbin/applications/MultiPhase/ex10/output/marker0.csv', '/home/eaulisa/FEMuS/femusbin/applications/MultiPhase/ex10/output/marker1.csv'])

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024

# show data in view
markerDisplay = Show(marker, spreadSheetView1, 'SpreadSheetRepresentation')

# get layout
layout1 = GetLayoutByName("Layout #1")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=0)

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(registrationName='TableToPoints1', Input=marker)
tableToPoints1.XColumn = 'X'
tableToPoints1.YColumn = 'Y'
tableToPoints1.ZColumn = 'Z'

# show data in view
tableToPoints1Display = Show(tableToPoints1, spreadSheetView1, 'SpreadSheetRepresentation')

# hide data in view
Hide(marker, spreadSheetView1)

# update the view to ensure updated data information
spreadSheetView1.Update()

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(tableToPoints1)

# show data in view
tableToPoints1Display_1 = Show(tableToPoints1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
tableToPoints1Display_1.Representation = 'Surface'
tableToPoints1Display_1.ColorArrayName = [None, '']
tableToPoints1Display_1.SelectTCoordArray = 'None'
tableToPoints1Display_1.SelectNormalArray = 'None'
tableToPoints1Display_1.SelectTangentArray = 'None'
tableToPoints1Display_1.OSPRayScaleArray = 'Nx'
tableToPoints1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
tableToPoints1Display_1.SelectOrientationVectors = 'None'
tableToPoints1Display_1.ScaleFactor = 0.0649
tableToPoints1Display_1.SelectScaleArray = 'None'
tableToPoints1Display_1.GlyphType = 'Arrow'
tableToPoints1Display_1.GlyphTableIndexArray = 'None'
tableToPoints1Display_1.GaussianRadius = 0.003245
tableToPoints1Display_1.SetScaleArray = ['POINTS', 'Nx']
tableToPoints1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
tableToPoints1Display_1.OpacityArray = ['POINTS', 'Nx']
tableToPoints1Display_1.OpacityTransferFunction = 'PiecewiseFunction'
tableToPoints1Display_1.DataAxesGrid = 'GridAxesRepresentation'
tableToPoints1Display_1.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tableToPoints1Display_1.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tableToPoints1Display_1.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# set active view
SetActiveView(spreadSheetView1)

# hide data in view
Hide(tableToPoints1, spreadSheetView1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# close an empty frame
layout1.Collapse(2)

# set active view
SetActiveView(renderView1)

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1330, 950)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [0.013399999999999995, 0.044999999999999984, 2.507577452241883]
renderView1.CameraFocalPoint = [0.013399999999999995, 0.044999999999999984, 0.0]
renderView1.CameraParallelScale = 0.6605555847584331

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).