#define CATCH_CONFIG_RUNNER


void runCatchTests(int argc, char *const argv[]) {
    Catch::Session().run(argc, argv);
}

int main(int argc, char *const argv[]) {
    runCatchTests(argc, argv);
    return 0;
}