package simulator

import "container/list"

const (
	// Key event callback type 'KeyDown' indicates an event that should be called
	// when a key is pressed down.
	KeyDown = iota

	// Key event callback type 'KeyUp' indicates an event that should be called
	// when a key is released.
	KeyUp

	// Key event callback type 'KeyHeld' indicates an event that should be called
	// before every single frame while a given key is being held down.
	KeyHeld
)

var (
	// Stores callbacks for possible key presses.
	callbacks = make(map[int]func(int))

	// Stores callbacks that must be called each frame. This map is updated
	// dynamically, so everything in it must be called every frame. When a callback
	// is no longer necessary, it gets deleted from the map.
	frameCallbacks = make(map[int]func())

	// A queue of callbacks that must be called on the next frame. When a frame is handled,
	// the queue is cleared and every callback is called and removed from the queue.
	callbackQueue = list.New()
)

/*** Exposed functions ***/

// Register a key event callback of a given type.
// Available types are:
//      KeyUp:    called when a key is released
//      KeyDown:  called when a key is pressed
//      KeyHeld:  called every frame while a key is down
func RegisterKey(key, eventType int, callback func()) {
	switch eventType {
	case KeyDown:
		registerKeyStateChanged(key, func(state int) {
			// Only respond to KeyDown events when a key is pressed down
			if state == 1 {
				callback()
			}
		})

	case KeyUp:
		registerKeyStateChanged(key, func(state int) {
			// Only respond to KeyUp events when a key is released
			if state == 0 {
				callback()
			}
		})

	case KeyHeld:
		registerKeyStateChanged(key, func(state int) {
			// When a key is pressed down, store the callback in the list of
			// callbacks that are called every frame.
			if state == 1 {
				frameCallbacks[key] = callback
			} else {
				// When a key is released, remove the callback from the list
				delete(frameCallbacks, key)
			}
		})
	}
}

/*** Internal functions ***/

// Register a callback to trigger when a key is pressed or unpressed.
// The callback receives the key state as an argument.
func registerKeyStateChanged(key int, callback func(int)) {
	callbacks[key] = callback
}

// Record input events from OpenGL key events. Call HandleInput() to cause
// event callbacks to happen. 
func onKey(key, state int) {
	// If we have a callback for this key, put it in the callback queue.
	if callback, ok := callbacks[key]; ok {
		callbackQueue.PushBack(func() {
			callback(state)
		})
	}
}

// Execute all event callbacks that have occurred since the last time
// HandleInput() was called and all the event callbacks that must be 
// called every frame for some duration.
func handleInput() {
	// Call all the callbacks that must be called every frame.
	for _, callback := range frameCallbacks {
		callback()
	}

	// Call the callbacks that are in the callback queue.
	for e := callbackQueue.Front(); e != nil; e = e.Next() {
		e.Value.(func())()
		callbackQueue.Remove(e)
	}
}
