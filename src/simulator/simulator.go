package simulator

import (
	"log"
	"math"
	"os"
	"time"
    "vector"

	"github.com/go-gl/gl"
	"github.com/go-gl/glfw"
	"github.com/go-gl/gltext"
	"github.com/go-gl/glu"
	"github.com/rhencke/glut"
)

const (
	// Motion rates
	moveRate   = 0.1
	rotateRate = 2.0

	// Camera positioning parameters
	initPhi      = 20.0
	initTheta    = 40.0
	initDistance = 6.0
	minDistance  = 1.0
	maxDistance  = 10
)

var (

	// Main loop will exit if running becomes false.
	running = true

	// Pause simulation but not visualization
	paused = false

	// Use help overlay
	helpOverlay = false

	// Display coordinate axes
	globalAxes = false

	// Enable wireframe drawing
	wireframe = false

	// Whether lighting is enabled
	lighting = true

	// Vertices and faces to draw
	meshes = make(map[string]*Mesh)

	// Textures to use
	// Index 0: ground texture
	// Index 1: wall texture
	textures = make([]gl.Texture, 2)

	// Whether to use texturing
	texturing = true

	// Initial translation
	phi      = initPhi
	theta    = initTheta
	distance = initDistance

	// Help strings; modifiable via AddHotkeyHelpText()
	keyboardHelpStrings = map[string]string{
		"ESC or Q": "Exit the simulation.",
		"W/S":      "Move closer or farther.",
		"A/D":      "Rotate around the simulation.",
		"R/F":      "Increase/decrease zenith angle.",
		"L":        "Toggle lighting.",
		"T":        "Toggle textures.",
		"V":        "Toggle wireframe models.",
		"H":        "Toggle the help menu.",
		"X":        "Toggle global coordinate axes.",
		"Space":    "Reset to original camera position.",
	}

	keyboardHelpLocs = map[string]int{
		"ESC or Q": 0,
		"W/S":      1,
		"A/D":      2,
		"R/F":      3,
		"L":        4,
		"T":        5,
		"V":        6,
		"H":        7,
		"X":        8,
		"Space":    9,
	}
)

var width, height, framesPerSecond int
var ticker <-chan time.Time

var font *gltext.Font
var boldFont *gltext.Font

/*** Exposed functions ***/

// Initialize the simulator, given the window name, width, height, and FPS.
func Init(name string, w, h, fps int) error {
	// Store simulation parameters
	width = w
	height = h
	framesPerSecond = fps

	// Compute frame duration from constant FPS, and initialize a timer.
	ticker = time.Tick(time.Duration(1000/framesPerSecond) * time.Millisecond)

	// Initialize GLFW.
	if err := glfw.Init(); err != nil {
		return err
	}

	// Open and initialize a non-resizable new window.
	glfw.OpenWindowHint(glfw.WindowNoResize, gl.TRUE)
	if err := glfw.OpenWindow(width, height, 8, 8, 8, 8, 0, 8, glfw.Windowed); err != nil {
		return err
	}
	glfw.SetSwapInterval(1)
	glfw.SetWindowTitle(name)

	// Load resources
	initResources()

	// Initialize OpenGL and all input action callbacks.
	initGL()
	initInput()

	return nil
}

// Exit the simulation, releasing resources nicely.
func Terminate() {
	glfw.CloseWindow()
	glfw.Terminate()

	// Release resources.
	font.Release()
	boldFont.Release()
}

// Returns true if the simulation is still running and being shown.
func Running() bool {
	return running && glfw.WindowParam(glfw.Opened) == 1
}

// Returns true if the simulation is paused
func Paused() bool {
	return paused
}

// Blocks until it is time for the next frame.
// Call this at the beginning of your main loop to set an FPS.
func WaitForNextFrame() {
	<-ticker
}

// Update the simulator; this should be called between WaitForNextFrame and Draw.
func Update() {
	// In each iteration of the loop, deal with input events since the last
	// frame, update the simulation, and redraw the window.
	handleInput()
}

// Update the rendering of the simulation.
func Draw() {
	// Clear the screen and depth buffer, reset matrices.
	gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)

	// Enable perspective for 3D drawing
	gl.MatrixMode(gl.PROJECTION)
	gl.LoadIdentity()
	glu.Perspective(45.0, float64(width)/float64(height), 0.1, 100.0)
	gl.MatrixMode(gl.MODELVIEW)
	gl.LoadIdentity()

	// Draw 3D objects
	drawScene()

	// Switch to orthographic mode for 2D overlays
	gl.MatrixMode(gl.PROJECTION)
	gl.LoadIdentity()
	gl.Ortho(0.0, float64(width), 0, float64(height), 0.0, 30.0)
	gl.MatrixMode(gl.MODELVIEW)
	gl.LoadIdentity()

	// Overlay 2D objects on top of the scene
	drawOverlays()

	// Draw buffer to screen.
	glfw.SwapBuffers()
}

// Add a line to the help text that you see when you press 'H'
func AddHotkeyHelpText(hotkey, description string) {
	keyboardHelpLocs[hotkey] = len(keyboardHelpLocs)
	keyboardHelpStrings[hotkey] = description
}

// Give textures to the ground and the walls of the simulation boundaries.
func TextureBoundaries(ground, walls string) {
	loadTexture(textures[0], ground)
	loadTexture(textures[1], walls)
}

// Add a mesh to be drawn. If this mesh has the same name as a mesh already in
// the simulation, the mesh will replace the older one.
func AddMesh(mesh *Mesh) {
	meshes[mesh.Name] = mesh
}

// Remove a mesh that is being drawn, given the mesh name.
func DeleteMesh(name string) {
	delete(meshes, name)
}

/*** Internal functions ***/

func initResources() {
	font = loadFont("data/UbuntuFont.ttf", 18)
	boldFont = loadFont("data/UbuntuBold.ttf", 22)

	// Allocate space for textures before loading them
	gl.GenTextures(textures)
}

func initInput() {
	// Exit the program with ESC or q.
	exit := func() { running = false }
	RegisterKey(glfw.KeyEsc, KeyDown, exit)
	RegisterKey('Q', KeyDown, exit)

	// Toggle help dialog with 'H'
	RegisterKey('H', KeyDown, func() {
		helpOverlay = !helpOverlay
	})

	// Toggle axes with 'X'
	RegisterKey('X', KeyDown, func() {
		globalAxes = !globalAxes
	})

	// Toggle whether lighting is enabled with l.
	RegisterKey('L', KeyDown, func() {
		if lighting = !lighting; lighting {
			gl.Enable(gl.LIGHTING)
		} else {
			gl.Disable(gl.LIGHTING)
		}
	})

	RegisterKey('T', KeyDown, func() {
		if texturing = !texturing; texturing {
			gl.Enable(gl.TEXTURE_2D)
		} else {
			gl.Disable(gl.TEXTURE_2D)
		}
	})

	// Toggle whether to draw models as wireframes.
	RegisterKey('V', KeyDown, func() {
		wireframe = !wireframe
	})

	// Movement controls
	RegisterKey('W', KeyHeld, func() {
		distance -= moveRate
		if distance < minDistance {
			distance = minDistance
		}
	})
	RegisterKey('S', KeyHeld, func() {
		distance += moveRate
		if distance > maxDistance {
			distance = maxDistance
		}
	})
	RegisterKey('R', KeyHeld, func() {
		theta -= rotateRate
		if theta < 0 {
			theta = 0
		}
	})
	RegisterKey('F', KeyHeld, func() {
		theta += rotateRate

		// Do not allow camera to be perfectly flat to avoid seeing below the ground
		if theta > 89.5 {
			theta = 89.5
		}
	})
	RegisterKey('D', KeyHeld, func() {
		phi += rotateRate
		if phi > 360 {
			phi -= 360
		}
	})
	RegisterKey('A', KeyHeld, func() {
		phi -= rotateRate
		if phi < 0 {
			phi += 360
		}
	})
	RegisterKey(' ', KeyHeld, func() {
		phi = initPhi
		theta = initTheta
		distance = initDistance
	})

	// Pause/unpause simulation
	pause := func() {
		paused = !paused
	}
	RegisterKey('p', KeyDown, pause)
	RegisterKey('P', KeyDown, pause)
	AddHotkeyHelpText("P", "Toggle simulation pause")


	// Start listening for key presses
	glfw.SetKeyCallback(onKey)
}

func initGL() {
	// Initialize the viewport and perspective.
	gl.Viewport(0, 0, width, height)

	gl.ShadeModel(gl.SMOOTH)
	gl.ClearColor(0, 0, 0, 0)
	gl.ClearDepth(1)
	gl.Enable(gl.DEPTH_TEST)
	gl.DepthFunc(gl.LEQUAL)
	gl.Hint(gl.PERSPECTIVE_CORRECTION_HINT, gl.NICEST)

	mat_specular := []float32{1.0, 1.0, 1.0, 1.0}
	mat_diffuse := []float32{0.5, 1.0, 0.5, 0.5}
	mat_shininess := []float32{50.0}
	light0_position := []float32{1.0, 1.0, 1.0, 0.0}
	light1_position := []float32{-1.0, -1.0, 1.0, 0.0}

	light1_diffuse := []float32{1.0, 1.0, 1.0, 1.0}
	light1_specular := []float32{0.1, 0.1, 0.1, 1.0}

	gl.Materialfv(gl.FRONT_AND_BACK, gl.SPECULAR, mat_specular)
	gl.Materialfv(gl.FRONT_AND_BACK, gl.SHININESS, mat_shininess)
	gl.Materialfv(gl.FRONT_AND_BACK, gl.DIFFUSE, mat_diffuse)
	gl.Lightfv(gl.LIGHT0, gl.POSITION, light0_position)
	gl.Lightfv(gl.LIGHT0, gl.SPECULAR, light1_specular)

	gl.Lightfv(gl.LIGHT1, gl.POSITION, light1_position)
	gl.Lightfv(gl.LIGHT1, gl.DIFFUSE, light1_diffuse)
	gl.Lightfv(gl.LIGHT1, gl.SPECULAR, light1_specular)

	gl.Enable(gl.LIGHTING)
	gl.Enable(gl.LIGHT0)
	gl.Enable(gl.LIGHT1)

	gl.Enable(gl.TEXTURE_2D)
}

func drawScene() {
	// Compute camera location based on spherical coordinates
	camera := vector.Vector{distance * math.Sin(Radian(theta)) * math.Cos(Radian(phi)),
		distance * math.Sin(Radian(theta)) * math.Sin(Radian(phi)),
		distance * math.Cos(Radian(theta))}

	// Compute "up" vector which defines camera orientation. We do this 
	// by taking the cross product of the camera vector with the unit vector
	// pointing in the direction of maximal increase of zenith.
	azimuthGrad := vector.Vector{-math.Sin(Radian(phi)), math.Cos(Radian(phi)), 0}
	up := camera.CrossProduct(azimuthGrad)

	// Set up camera direction
	gl.MatrixMode(gl.MODELVIEW)
	gl.LoadIdentity()
	glu.LookAt(camera.X, camera.Y, camera.Z, 0, 0, 0.5, up.X, up.Y, up.Z)

	// Disable lighting when using wireframing
	if lighting && wireframe {
		gl.Disable(gl.LIGHTING)
		defer gl.Enable(gl.LIGHTING)
	}

	// Draw simulation boundaries to avoid ugly black background
	drawGround()

	// Draw all meshes
	for _, mesh := range meshes {
		drawMesh(mesh)
	}

	if globalAxes {
		drawCoordinateAxes()
	}
}

func drawMesh(mesh *Mesh) {
	vertices := mesh.Vertices

	// Apply local transformations
	// TODO: Fix rotation of meshes. Need to convert from Euler angles to axis-angle to do rotation in GL.
	gl.PushMatrix()
	gl.Translatef(mesh.dx, mesh.dy, mesh.dz)
	gl.Rotatef(mesh.rotZ, 0, 0, 1)
	gl.Rotatef(mesh.rotY, 0, 1, 0)
	gl.Rotatef(mesh.rotX, 1, 0, 0)
	gl.Scalef(mesh.sx, mesh.sy, mesh.sz)
	defer gl.PopMatrix()

	// Draw each face separately.
	// TODO: Fix this to use vertex arrays.
	for _, face := range mesh.Faces {
		if wireframe {
			gl.Color3f(1.0, 1.0, 1.0)
			gl.Begin(gl.LINE_LOOP)
		} else if len(face) == 4 {
			gl.Begin(gl.QUADS)
		} else if len(face) == 3 {
			gl.Begin(gl.TRIANGLES)
		} else {
			panic("Unknown GL polygon type.")
		}

		for _, vertIndex := range face {
			vertex := vertices[vertIndex]
			gl.Normal3f(float32(mesh.VertexNormals[vertIndex].X), float32(mesh.VertexNormals[vertIndex].Y), float32(mesh.VertexNormals[vertIndex].Z))
			gl.Vertex3f(float32(vertex.X), float32(vertex.Y), float32(vertex.Z))
		}

		gl.End()
	}
}

func drawOverlays() {
	// Disable lighting for overlays
	if lighting {
		gl.Disable(gl.LIGHTING)
		defer gl.Enable(gl.LIGHTING)
	}

	if helpOverlay {
		overlayBorder := float32(30)

		// Allow blending for semi-transparent overview
		gl.Enable(gl.BLEND)
		defer gl.Disable(gl.BLEND)

		// Draw transparent quad covering most of scene
		gl.Color4f(0.1, 0.1, 0.1, 0.9)
		gl.Begin(gl.QUADS)
		gl.Vertex2f(overlayBorder, overlayBorder)
		gl.Vertex2f(float32(width)-overlayBorder, overlayBorder)
		gl.Vertex2f(float32(width)-overlayBorder, float32(height)-overlayBorder)
		gl.Vertex2f(overlayBorder, float32(height)-overlayBorder)
		gl.End()

		// Text color
		gl.Color3f(0.9, 0.9, 0.9)

		// Help text heading
		heading := "Keyboard Commands"
		headingWidth, _ := boldFont.Metrics(heading)
		boldFont.Printf(float32(width-headingWidth)/2, float32(overlayBorder+10), heading)

		// Overlay help text
		borderPadding := float32(50)
		lineSpacing := float32(30)
		tabLocation := float32(300)
		for key, str := range keyboardHelpStrings {
			i := keyboardHelpLocs[key]
			x := overlayBorder + borderPadding
			y := overlayBorder + borderPadding + lineSpacing*float32(i)
			font.Printf(float32(x), float32(y), "%s", key)
			font.Printf(tabLocation, float32(y), "%s", str)
		}
	}
}

func drawCoordinateAxes() {
	// Draw coordinate axes and arrows in colors, without lighting.
	if lighting {
		gl.Disable(gl.LIGHTING)
		defer gl.Enable(gl.LIGHTING)
	}

	// Red x-axis line.
	gl.Color3f(1.0, 0.0, 0.0)
	gl.Begin(gl.LINES)
	gl.Vertex3f(0.0, 0.0, 0.0)
	gl.Vertex3f(1.0, 0.0, 0.0)
	gl.End()

	// Red arrow.
	gl.PushMatrix()
	gl.Translatef(1.0, 0.0, 0.0)
	gl.Rotatef(90, 0.0, 1.0, 0.0)
	glut.SolidCone(0.04, 0.2, 10, 10)
	gl.PopMatrix()

	// Green y-axis line.
	gl.Color3f(0.0, 1.0, 0.0)
	gl.Begin(gl.LINES)
	gl.Vertex3f(0.0, 0.0, 0.0)
	gl.Vertex3f(0.0, 1.0, 0.0)
	gl.End()

	// Green arrow.
	gl.PushMatrix()
	gl.Translatef(0.0, 1.0, 0.0)
	gl.Rotatef(-90, 1.0, 0.0, 0.0)
	glut.SolidCone(0.04, 0.2, 10, 10)
	gl.PopMatrix()

	// Blue z-axis line.
	gl.Color3f(0.0, 0.0, 1.0)
	gl.Begin(gl.LINES)
	gl.Vertex3f(0.0, 0.0, 0.0)
	gl.Vertex3f(0.0, 0.0, 1.0)
	gl.End()

	// Blue arrow.
	gl.PushMatrix()
	gl.Translatef(0.0, 0.0, 1.0)
	gl.Rotatef(-90, 0.0, 0.0, 1.0)
	glut.SolidCone(0.04, 0.2, 10, 10)
	gl.PopMatrix()
}

func drawGround() {
	// Use a white material to use the texture's natural color.
	mat_specular := []float32{
		1.0, 1.0, 1.0, 1.0}
	mat_diffuse := []float32{
		1.0, 1.0, 0.8, 0.5}
	mat_shininess := []float32{
		100.0}
	gl.Materialfv(gl.FRONT_AND_BACK, gl.SPECULAR, mat_specular)
	gl.Materialfv(gl.FRONT_AND_BACK, gl.SHININESS, mat_shininess)
	gl.Materialfv(gl.FRONT_AND_BACK, gl.DIFFUSE, mat_diffuse)

	// Draw ground quad
	textures[0].Bind(gl.TEXTURE_2D)
	roomSize := float32(maxDistance * 1.1)
	gl.Begin(gl.QUADS)
	gl.Normal3f(0, 0, 1)
	gl.TexCoord2f(0, 0)
	gl.Vertex3f(-roomSize, -roomSize, 0)
	gl.TexCoord2f(0, 1)
	gl.Vertex3f(-roomSize, roomSize, 0)
	gl.TexCoord2f(1, 1)
	gl.Vertex3f(roomSize, roomSize, 0)
	gl.TexCoord2f(1, 0)
	gl.Vertex3f(roomSize, -roomSize, 0)
	gl.End()
	textures[0].Unbind(gl.TEXTURE_2D)

	// Draw walls
	textures[1].Bind(gl.TEXTURE_2D)
	gl.Begin(gl.QUADS)
	gl.Normal3f(1, 0, 0)
	gl.TexCoord2f(0, 0)
	gl.Vertex3f(-roomSize, -roomSize, 0)
	gl.TexCoord2f(0, 1)
	gl.Vertex3f(-roomSize, roomSize, 0)
	gl.TexCoord2f(1, 1)
	gl.Vertex3f(-roomSize, roomSize, roomSize)
	gl.TexCoord2f(1, 0)
	gl.Vertex3f(-roomSize, -roomSize, roomSize)

	gl.Normal3f(-1, 0, 0)
	gl.TexCoord2f(0, 0)
	gl.Vertex3f(roomSize, -roomSize, 0)
	gl.TexCoord2f(0, 1)
	gl.Vertex3f(roomSize, roomSize, 0)
	gl.TexCoord2f(1, 1)
	gl.Vertex3f(roomSize, roomSize, roomSize)
	gl.TexCoord2f(1, 0)
	gl.Vertex3f(roomSize, -roomSize, roomSize)

	gl.Normal3f(0, 1, 0)
	gl.TexCoord2f(0, 0)
	gl.Vertex3f(-roomSize, -roomSize, 0)
	gl.TexCoord2f(0, 1)
	gl.Vertex3f(roomSize, -roomSize, 0)
	gl.TexCoord2f(1, 1)
	gl.Vertex3f(roomSize, -roomSize, roomSize)
	gl.TexCoord2f(1, 0)
	gl.Vertex3f(-roomSize, -roomSize, roomSize)

	gl.Normal3f(0, -1, 0)
	gl.TexCoord2f(0, 0)
	gl.Vertex3f(-roomSize, roomSize, 0)
	gl.TexCoord2f(0, 1)
	gl.Vertex3f(roomSize, roomSize, 0)
	gl.TexCoord2f(1, 1)
	gl.Vertex3f(roomSize, roomSize, roomSize)
	gl.TexCoord2f(1, 0)
	gl.Vertex3f(-roomSize, roomSize, roomSize)
	gl.End()

	textures[1].Unbind(gl.TEXTURE_2D)
}

func loadFont(filename string, size int) *gltext.Font {
	fd, err := os.Open(filename)
	if err != nil {
		log.Printf("Error opening file %s: %v", filename, err)
		os.Exit(1)
	}
	defer fd.Close()

	font, err := gltext.LoadTruetype(fd, int32(size), 32, 127, gltext.LeftToRight)
	if err != nil {
		log.Printf("Error loading font: %v", err)
		os.Exit(1)
	}

	return font
}

func loadTexture(texture gl.Texture, filename string) {
	// Load the ground texture. Note that LoadTexture2D loads the image data
	// from a TGA file into whatever texture is bound at the moment.
	texture.Bind(gl.TEXTURE_2D)
	if !glfw.LoadTexture2D(filename, 0) {
		log.Printf("Failed to load texture %s.", filename)
		os.Exit(1)
	}

	// Set the texture minifaction and magnification filters. These specify how
	// to get a value at a given point when the level of detail in the texture
	// is too low or too high. gl.LINEAR uses linear interpolation between
	// nearby texture points; gl.NEAREST just uses the nearest texture point.
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR)
}
