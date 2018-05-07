#ifndef SEC
#define SEC

static Curve E_0 = {    { {0x5CF98103A5247EFD, 0x7}, {0xE8C77E323E882E70, 0xF} 
},
                        { {0x8A2E5D3757296583, 0x2}, {0x209398683CC050BE, 0x10} 
}
                   };
static Point P_a = {    { {0x27DCDDAFF812E2F4, 0x9}, {0xE28C6D8C26D05B9D, 0x8} 
},
                        { {0xDE8598E7E4633A84, 0xB}, {0x771DA880B2FC22D, 0x10} 
},
                        { {0x1, 0x0}, {0x0, 0x0} }
                   };
static Point Q_a = {    { {0xCB863062A3934C5C, 0x9}, {0xE2C400C5D2EC6CB2, 0xA} 
},
                        { {0x22173C0938306C89, 0x3}, {0x7EA4F807C9D80010, 0x1} 
},
                        { {0x1, 0x0}, {0x0, 0x0} }
                   };
static Point P_b = {    { {0x8495C61DC1B71B6, 0x0}, {0x7971C69C8547B66B, 0xC} },
                        { {0xE284471C9F5E3BB0, 0x0}, {0x9916A5727AA092BF, 0xE} 
},
                        { {0x1, 0x0}, {0x0, 0x0} }
                   };
static Point Q_b = {    { {0x8A49E0C68B0AF958, 0x6}, {0xD951EA62F7387E23, 0xB} 
},
                        { {0xC4C370348226B79, 0x6}, {0x5E20B9FAE60BD694, 0xB} },
                        { {0x1, 0x0}, {0x0, 0x0} }
                   };
//---------------------------------------------
static Curve E_A = {    { {0xE0E92F431EE9F99, 0xC}, {0x5EF60ECBA45DFFA7, 0x6} },
                        { {0x26629EE96A72083, 0xB}, {0x9CEC01B5D21F0BFF, 0x6} }
                   };
static Point phi_P_b = { { {0x8FF89139011713D9, 0xC}, {0x7C5CB28DEF61494A, 0x9} 
},
                         { {0xD1FC32D815775342, 0xB}, {0xBA0522F6628FDE08, 0x4} 
},
                        { {0x1, 0x0}, {0x0, 0x0} }
                   };
static Point phi_Q_b = { { {0xCD6EB12768A418D0, 0x8}, {0x36F073DD956C1169, 0x4} 
},
                         { {0x7A7413B77977EBB1, 0x7}, {0x9588E881DAEC0CC6, 0x6} 
},
                        { {0x1, 0x0}, {0x0, 0x0} }
                   };
static Point S_a = {    { {0xE3DD4987C45C8306, 0x0}, {0x8D10D95C66B339C, 0x9} },
                        { {0x59D13209B9F1FCCC, 0x4}, {0x17D0F40238A60A22, 0x3} 
},
                        { {0x1, 0x0}, {0x0, 0x0} }
                   };
static Point T_a = {    { {0x525FC6B4EF8B1992, 0xD}, {0xD5E896D79ED5924C, 0x4} 
},
                        { {0x14D4B381DEA3CF39, 0xC}, {0x67E6B3A91786C6DC, 0x10} 
},
                        { {0x1, 0x0}, {0x0, 0x0} }
                   };
//---------------------------------------------
static Curve E_B = {    { {0xEC81B92C9E10F766, 0x5}, {0xAB26F7D9275D0ABB, 0x6} 
},
                        { {0xEEDD9DC3AFD0F90B, 0x5}, {0x31FF77F38857AB17, 0x12} 
}
                   };
static Point psi_P_a = { { {0xF7A68203E18681CD, 0x6}, {0x4DFFA0FBA84C5371, 0x1} 
},
                         { {0xABFBA23240304A5, 0x7}, {0xCDEB2D61E5D62C0A, 0xA} 
},
                        { {0x1, 0x0}, {0x0, 0x0} }
                   };
static Point psi_Q_a = { { {0xA64A5BD16C24813E, 0xA}, {0x2FBD08E80A931C57, 0xA} 
},
                         { {0x9B3CE7F07687AA3A, 0x6}, {0xE7C127B394BE4D7E, 0xC} 
},
                        { {0x1, 0x0}, {0x0, 0x0} }
                   };
static Point S_b = {    { {0xCCE21C3AD3D06A5A, 0x4}, {0xD45F3C7AFCF9EA91, 0xF} 
},
                        { {0x725513C5A0609D71, 0x3}, {0xD65F01B869941BFE, 0x6} 
},
                        { {0x1, 0x0}, {0x0, 0x0} }
                   };
static Point T_b = {    { {0x4469FCE20DC1F079, 0xE}, {0x2AF00671912C4BF3, 0x3} 
},
                        { {0xBCBF753EA30C8B38, 0x2}, {0x25BB2D5B82E57C9B, 0x9} 
},
                        { {0x1, 0x0}, {0x0, 0x0} }
                   };
#endif
