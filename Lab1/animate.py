from celluloid import Camera
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from ordered_set import OrderedSet


def parse_file(file_name: str = "output", every_n: int = 1) -> dict[str: list[tuple[float, float]]]:
    bodies_trajectories = {}
    with open(file_name, "r", encoding='utf-8') as file:
        rows = [line_stripped.lower() for line in file.readlines()[1:] if (line_stripped := line.strip())]

        for row in rows:
            if "cycle" in row:
                continue
            tokens = [token.strip() for token in row.split()]
            body_key = row.split(" :")[0]
            bodies_trajectories.setdefault(body_key, OrderedSet()).add((float(tokens[3]), float(tokens[4])))

        for body, trajectory in bodies_trajectories.items():
            bodies_trajectories[body] = [trajectory[i] for i in range(0, len(trajectory), every_n)]

    return bodies_trajectories


def normolize_points(bodies_trajectories) -> dict[str: list[tuple[float, float]]]:
    bodies_trajectories_normolized = {}
    for body, trajectory in bodies_trajectories.items():
        x_points = [point[0] for point in trajectory]
        y_points = [point[1] for point in trajectory]
        norm_x = [float(i)/sum(x_points) for i in x_points]
        norm_y = [float(i)/sum(y_points) for i in y_points]

        norm_trajectory = [(x, y) for x, y in zip(norm_x, norm_y)]

        bodies_trajectories_normolized[body] = norm_trajectory
    return bodies_trajectories_normolized


def generate_animation(bodies_trajectories: dict, animation_file: str = "animation.gif", max_points=1000):
    fig, ax = plt.subplots()
    camera = Camera(fig)
    points_n = len(list(bodies_trajectories.values())[0])
    # bodies = bodies_trajectories.keys()
    for point_i in range(points_n):
        for (body, trajectory), color in zip(bodies_trajectories.items(), mcolors.TABLEAU_COLORS.values()):
            # plt.plot(trajectory[point_i][0], trajectory[point_i][1], 'bo', color=color, lw=0.8)
            # ax.annotate(txt, (z[i], y[i]))
            # x = [trajectory[point_i][0] for trajectory]
            ax.scatter(trajectory[point_i][0], trajectory[point_i][1], color=color)
            ax.annotate(body, (trajectory[point_i][0], trajectory[point_i][1]), fontsize=7)

        camera.snap()

    animation = camera.animate(repeat=True)

    animation.save(animation_file, )
    print(f"Animation saved to {animation_file}")


def main():
    bodies_trajectories = parse_file(file_name="obj_16", every_n=500)
    # bodies_trajectories = normolize_points(bodies_trajectories)
    # print(bodies_trajectories.keys())
    generate_animation(bodies_trajectories)


if __name__ == '__main__':
    main()
