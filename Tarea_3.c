{
  "nbformat": 4,
  "nbformat_minor": 0,
  "metadata": {
    "colab": {
      "provenance": [],
      "authorship_tag": "ABX9TyNXHaSjV7dMJIEHstT0IUlu",
      "include_colab_link": true
    },
    "kernelspec": {
      "name": "python3",
      "display_name": "Python 3"
    },
    "language_info": {
      "name": "python"
    }
  },
  "cells": [
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "view-in-github",
        "colab_type": "text"
      },
      "source": [
        "<a href=\"https://colab.research.google.com/github/DCI-alxogm/me2025-clase-rocillorn2021-create/blob/main/Tarea_3.c\" target=\"_parent\"><img src=\"https://colab.research.google.com/assets/colab-badge.svg\" alt=\"Open In Colab\"/></a>"
      ]
    },
    {
      "cell_type": "code",
      "execution_count": 2,
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "xh_tdDb1L2ld",
        "outputId": "95872b84-483f-4808-9b19-ab40e957f981"
      },
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "/usr/bin/ld: /tmp/ccPfmKn8.o: in function `main':\n",
            "Tarea_3.c:(.text+0x1c8): undefined reference to `pow'\n",
            "/usr/bin/ld: Tarea_3.c:(.text+0x287): undefined reference to `pow'\n",
            "/usr/bin/ld: Tarea_3.c:(.text+0x373): undefined reference to `pow'\n",
            "collect2: error: ld returned 1 exit status\n"
          ]
        }
      ],
      "source": [
        "!gcc Tarea_3.c -o -lm Tarea_3.o"
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "!gcc Tarea_3.c -o -lm Tarea_3.o"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "KZOKjZDnNgUp",
        "outputId": "22e84453-e725-4a8e-ae1f-a0e6c7020eb2"
      },
      "execution_count": 3,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stdout",
          "text": [
            "/usr/bin/ld: /tmp/ccFBHsDr.o: in function `main':\n",
            "Tarea_3.c:(.text+0x1c8): undefined reference to `pow'\n",
            "/usr/bin/ld: Tarea_3.c:(.text+0x287): undefined reference to `pow'\n",
            "/usr/bin/ld: Tarea_3.c:(.text+0x373): undefined reference to `pow'\n",
            "collect2: error: ld returned 1 exit status\n"
          ]
        }
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "!gcc Tarea_3.c -o Tarea_3.o -lm\n"
      ],
      "metadata": {
        "id": "lgfmE4ohOF54"
      },
      "execution_count": 4,
      "outputs": []
    }
  ]
}