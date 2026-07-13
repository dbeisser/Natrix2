#!/bin/bash

# Build and upload the Natrix2 Docker image.
# Only for the Natrix2 dev branch.

set -u

IMAGE="dbeisser/natrix2:latest"
REQUIRED_BRANCH="dev"
REQUIRED_REMOTE="dbeisser/Natrix2"

# Locate project root
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(dirname "$script_dir")"

cd "$project_root" || {
    echo "ERROR: Cannot access the Natrix2 project root."
    exit 1
}

# Launcher
echo ""
echo "Natrix2 Docker Builder"
echo "Build and upload the Docker image"
echo "Docker image: $IMAGE"
echo "Required branch: $REQUIRED_BRANCH"

# Check dependencies
for command in git docker; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "ERROR: Required command '$command' is not installed."
        exit 1
    fi
done

# Verify repository
remote_url="$(git remote get-url origin 2>/dev/null || true)"

if [[ "$remote_url" != *"$REQUIRED_REMOTE"* ]]; then
    echo "ERROR: This script may only be used in the official Natrix2 repository."
    echo "Detected remote: ${remote_url:-none}"
    exit 1
fi

# Verify branch
current_branch="$(git branch --show-current 2>/dev/null)"

if [[ -z "$current_branch" ]]; then
    echo "ERROR: Cannot determine the current Git branch."
    exit 1
fi

if [[ "$current_branch" != "$REQUIRED_BRANCH" ]]; then
    echo "ERROR: Docker images may only be built from the '$REQUIRED_BRANCH' branch."
    echo "Current branch: $current_branch"
    echo "Switch branch with: git switch $REQUIRED_BRANCH"
    exit 1
fi

# Check Dockerfile
if [[ ! -f "Dockerfile" ]]; then
    echo "ERROR: Dockerfile not found in the Natrix2 project root."
    exit 1
fi

# Check Docker
if ! docker info >/dev/null 2>&1; then
    echo "ERROR: Docker is not running or cannot be accessed."
    exit 1
fi

# Read Docker Hub username
echo ""
echo "Docker Hub username:"
read -r -p "> " docker_username

if [[ -z "$docker_username" ]]; then
    echo "ERROR: No Docker Hub username entered."
    exit 1
fi

# Check for existing local image
SKIP_BUILD=false

if docker image inspect "$IMAGE" >/dev/null 2>&1; then
    echo ""
    echo "A local Docker image was found."
    read -r -p "Do you want to rebuild it? [yes/no]: " rebuild

    case "${rebuild,,}" in
        yes|y)
            ;;
        no|n|"")
            echo "Using existing Docker image."
            SKIP_BUILD=true
            ;;
        *)
            echo "ERROR: Invalid response."
            exit 1
            ;;
    esac
fi

# Build image if required
if [[ "$SKIP_BUILD" != true ]]; then
    echo ""
    echo "Building Docker image..."
    echo "Source branch: $current_branch"
    echo "Docker image:  $IMAGE"
    echo ""

    if ! docker build --pull --tag "$IMAGE" .; then
        echo ""
        echo "ERROR: Docker image build failed."
        exit 1
    fi
fi

# Confirm upload
echo ""
echo "Docker image built successfully."
echo "Docker image: $IMAGE"
echo ""

read -r -p "Upload '$IMAGE' to Docker Hub? [yes/no]: " confirmation

case "${confirmation,,}" in
    yes|y)
        ;;
    no|n|"")
        echo "Upload cancelled. The image remains available locally."
        exit 0
        ;;
    *)
        echo "ERROR: Invalid response. Upload cancelled."
        exit 1
        ;;
esac

# Log in to Docker Hub if required
echo ""
docker_config="${DOCKER_CONFIG:-$HOME/.docker}/config.json"

if [[ -f "$docker_config" ]] && grep -Eq '"(auths|credsStore|credHelpers)"' "$docker_config"; then
    echo "Docker Hub login detected."
else
    echo "No Docker Hub login detected."
    echo "Log in to Docker Hub as '$docker_username'."
    echo "Use your personal access token as the password."

    if ! docker login --username "$docker_username"; then
        echo "ERROR: Docker Hub login failed."
        exit 1
    fi
fi

# Upload image
echo ""
echo "Uploading Docker image..."

if ! docker push "$IMAGE"; then
    echo ""
    echo "ERROR: Docker image upload failed."
    echo "Verify that '$docker_username' has write access to '$IMAGE'."
    exit 1
fi

# Success
echo ""
echo "Docker image uploaded successfully."
echo "Docker image: $IMAGE"
