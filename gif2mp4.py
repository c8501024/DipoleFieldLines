import os
from moviepy import VideoFileClip


def convert_gif_to_mp4(gif_path):
    mp4_path = gif_path.rsplit(".", 1)[0] + ".mp4"
    print(f"Converting: {gif_path} -> {mp4_path}")
    clip = VideoFileClip(gif_path)
    clip.write_videofile(
        mp4_path,
        fps=clip.fps,  # FPS für die Konvertierung
        codec="libx264",
        audio=False,
        bitrate="1000k",
        ffmpeg_params=["-crf", "20", "-pix_fmt", "yuv420p"],
    )


def convert_all_gifs_in_dir():
    current_dir = os.getcwd()
    gif_files = [f for f in os.listdir(current_dir) if f.lower().endswith(".gif")]

    if not gif_files:
        print("No .gif files found in the current directory.")
        return

    for gif in gif_files:
        try:
            convert_gif_to_mp4(gif)
        except Exception as e:
            print(f"Failed to convert {gif}: {e}")


if __name__ == "__main__":
    convert_all_gifs_in_dir()
