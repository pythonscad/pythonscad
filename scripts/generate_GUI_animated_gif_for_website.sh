#!/bin/sh

# ==============================================================================
# ANIMATED GIF GENERATOR WITH SMOOTH MORPH TRANSITIONS
# ==============================================================================
# USAGE:
#   1. Place your source images in the same folder as this script.
#   2. Name them sequentially (e.g., frame1.png, frame2.png ... frame10.png).
#      Using two digits (frame01.png, frame02.png) ensures proper sorting order.
#   3. Make the script executable: chmod +x generate.sh
#   4. Run it: ./generate.sh
#
# CONFIGURATION:
#   Adjust the millisecond variables below to experiment with timing.
# ==============================================================================

STATIC_MS=1500
TRANSITION_TOTAL_MS=100
NUM_TRANSITION_FRAMES=3

# Calculate ticks (1 tick = 10ms for ImageMagick)
DELAY_STATIC=$(( STATIC_MS / 10 ))
DELAY_TRANS=$(( (TRANSITION_TOTAL_MS / NUM_TRANSITION_FRAMES) / 10 ))

# 1. Get a sorted list of all your frames
FILES=$(ls frame*.png | sort -V)
FIRST_FILE=$(echo "$FILES" | head -n 1)

# 2. Generate the morphed sequence
# FIX: We explicitly set '-delay $DELAY_TRANS' before morphing so the
# generated frames inherit the correct timing metadata.
convert -delay $DELAY_TRANS $FILES $FIRST_FILE -morph $NUM_TRANSITION_FRAMES morphed_%03d.png

# 3. Build the GIF command dynamically
cmd="convert"
count=0
step=$(( NUM_TRANSITION_FRAMES + 1 ))

# Loop through the generated morph frames
for file in morphed_*.png; do
    # Every 'step' frame is an original keyframe stage.
    # We explicitly override its delay to the static timing.
    if [ $(( count % step )) -eq 0 ]; then
        cmd="$cmd -delay $DELAY_STATIC $file"
    else
        # For intermediate frames, we explicitly pass the delay again to ensure override
        cmd="$cmd -delay $DELAY_TRANS $file"
    fi
    count=$(( count + 1 ))
done

# 4. Finalize, remove the very last duplicated keyframe, and compile.
cmd=$(echo "$cmd" | sed "s/-delay $DELAY_STATIC $file//")
$cmd -loop 0 animated_stages.gif

# 5. Clean up temporary files
rm morphed_*.png

echo "Done! Generated seamless loop in animated_stages.gif"
