TGraph* get_predictions_cgc() {
    TGraph* gCGC = new TGraph(21);
    gCGC->SetName("gCGC");
    gCGC->SetPoint(0, 3, 25.663);
    gCGC->SetPoint(1, 2.7, 25.397);
    gCGC->SetPoint(2, 2.4, 25.018);
    gCGC->SetPoint(3, 2.1, 24.546);
    gCGC->SetPoint(4, 1.8, 23.967);
    gCGC->SetPoint(5, 1.5, 23.276);
    gCGC->SetPoint(6, 1.2, 22.496);
    gCGC->SetPoint(7, 0.9, 21.692);
    gCGC->SetPoint(8, 0.6, 20.962);
    gCGC->SetPoint(9, 0.3, 20.456);
    gCGC->SetPoint(10, 0, 20.277);
    gCGC->SetPoint(11, -0.3, 20.362);
    gCGC->SetPoint(12, -0.6, 20.493);
    gCGC->SetPoint(13, -0.9, 20.454);
    gCGC->SetPoint(14, -1.2, 20.117);
    gCGC->SetPoint(15, -1.5, 19.546);
    gCGC->SetPoint(16, -1.8, 18.862);
    gCGC->SetPoint(17, -2.1, 18.294);
    gCGC->SetPoint(18, -2.4, 17.828);
    gCGC->SetPoint(19, -2.7, 17.345);
    gCGC->SetPoint(20, -3, 16.848);

    return gCGC;
}
