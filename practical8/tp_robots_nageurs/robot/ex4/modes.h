#ifndef __MODES_H
#define __MODES_H

/// Idle mode: do nothing
#define IMODE_IDLE          0

/// Motor move mode
#define IMODE_MOTOR_DEMO    1

/// The main loop for mode switching
void main_mode_loop(void);
void motor_demo_mode(void);

#endif // __MODES_H
